%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% Define constants
c=2.998E-5;% Speed of light cm/fs
q_e=1.60218e-19; % Fundamental charge of an electron in C
Dconv=3.33564e-30; % Conversion factor to D
kT=200; % kT in cm^-1

%% Default material parameters
% Solvent Raman differential cross section prefactor in m^2
K=1.92e-30; % Dichloromethane
% Solvent density in g/m^3
rho_Solv=1326600; % Dichloromethane
% Solvent molar mass in g/mol
mm_Solv=84.93; % Dichloromethane
% Solvent vibrational mode in cm^-1
wsolv=702; % Dichloromethane
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% refractive index Solute
n_w_t_Sol=1.39; % also just a number
% refractive index Solvent
n_w_t_Solv=1.424; % Dichloromethane
% Electronic energy gap origin Solute in cm^-1
weg=22720;
% Electronic energy gap origin Solvent in cm^-1
weg_Solv=126900; % Dichloromethane
% Signal frequency
w_t=weg;
% PNA parameters
[~,~,kappa,lambda,weg,translen,n_w_t,...
    wvibs,Disp]=PNA_parameters('dichloromethane');
% Electronic dephasing in cm^-1
gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
% Transition dipole in D
mu_eg=translen*1e-10*q_e/Dconv;

FSRS_Methanol_dts=10;

% Set number of vibrational modes
nmode=2;length(wvibs);
% Set number of vibrational quanta to apportion
nquanta=3;
% Compute Solvent concentration in mol m^-3
C_Solv=rho_Solv/mm_Solv;
C_Sol=2e-4;

%% Default laser paramters
% Spectral width actinic pump
LAMBDA_ap=210;
% Spectral width raman pump
LAMBDA_rp=40;
%
FSRS_Methanol_nt=128;

dt=FSRS_Methanol_dts;
w_ap=weg;
w_rp=w_ap;
wvib=0;
for ii=1:10
    Disp=Disp./Disp*ii*2.5/10;
    % Compute E(3):E(5) prefactor
    alpha2=polarizability(K,weg_Solv,w_ap);
    prefactor_OffRes=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_ap-wvib)/4/pi^2;
    prefactor_Res=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_ap-wvib);
    
    %% Compute matrix of basis states, vibrational energies, and overlap integrals
    % Basis states and vibrational energies
    [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
    % Overlap integrals
    [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);
    
    %% Consolidate parameters
    % Consolidate material parameters into structure variable
    % parameters_material
    parameters_material.disp=Disp;
    parameters_material.gamma_eg=gamma_eg;
    parameters_material.gamma_vib=gamma_vib;
    parameters_material.weg=weg;
    parameters_material.wvib=wvib;
    parameters_material.wsolv=wsolv;
    
    % Consolidate laser parameters into structure variable parameters_laser
    parameters_laser.LAMBDA_ap=LAMBDA_ap;
    parameters_laser.LAMBDA_rp=LAMBDA_rp;
    parameters_laser.w_ap=w_ap;
    parameters_laser.w_rp=w_rp;
    parameters_laser.dt=dt;
    parameters_laser.nt=FSRS_Methanol_nt;
    tic
    %% Compute third order cascade response functions
    [FSRS_Methanol_OffRes_Cascade,FSRS_Methanol_OffRes_Direct] = cascade_FSRS_OffResSpec_End(...
        wviball_TC,nmode*nquanta,ovlp_TC,parameters_laser,parameters_material);
    FSRS_Methanol_OffRes_Cascade=prefactor_OffRes*3e10*FSRS_Methanol_OffRes_Cascade;
    
    [FSRS_Methanol_Res_Cascade,FSRS_Methanol_Res_Direct(ii,:)] = cascade_FSRS_ResSpec_End(...
        wviball_TC,nmode*nquanta,ovlp_TC,parameters_laser,parameters_material);
    FSRS_Methanol_Res_Cascade=prefactor_Res/3e10*FSRS_Methanol_Res_Cascade;
    toc
    
end
figure;
hold on
mx=max(abs(real(FSRS_Methanol_Res_Direct)),[],2);
mxx=max(mx);
mx=abs(log10(mx/mxx))/3.2476;
for ii=1:10
    plot([0:1/c/128/10:127/128/c/10],...
        abs(real(FSRS_Methanol_Res_Direct)),'Color',[ii/10,0,0]);xlim([0 1600]);
end
% save('FSRS_Methanol','FSRS_Methanol_OffRes_Cascade','FSRS_Methanol_OffRes_Direct',...
%     'FSRS_Methanol_Res_Cascade','FSRS_Methanol_Res_Direct',...
%     'FSRS_Methanol_dts','FSRS_Methanol_nt');
'done'