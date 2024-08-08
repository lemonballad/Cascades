%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
% Speed of light cm/fs
c=2.998E-5;
q_e=1.60218e-19;
Dconv=3.33564e-30;
kT=200;

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
    wvib,Disp]=PNA_parameters('dichloromethane');
% Vibrational modes to be examined in cm^-1
wvibs=wvib([1 3]);
% Unitless mode displacements
disp_pna=Disp([1 3]);
% Electronic dephasing in cm^-1
gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
% Transition dipole in D
mu_eg=translen*1e-10*q_e/Dconv;
% Time step for specific vibrational mode.
Full_2dRR_dts=10;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;

% Compute Solvent concentration in mol m^-3
C_Solv=rho_Solv/mm_Solv;
C_Sol=2e-4;

Full_2dRR_nt=128;
Full_2dRR_OffRes_Cascade=complex(zeros(Full_2dRR_nt,Full_2dRR_nt,2,'double'));
Full_2dRR_OffRes_Direct=complex(zeros(Full_2dRR_nt,Full_2dRR_nt,2,'double'));
Full_2dRR_Res_Cascade=complex(zeros(Full_2dRR_nt,Full_2dRR_nt,2,'double'));
Full_2dRR_Res_Direct=complex(zeros(Full_2dRR_nt,Full_2dRR_nt,2,'double'));
matlabpool open 10
for iv=1:2
    %% Default laser paramters
    wvib=wvibs(iv);
    dt=Full_2dRR_dts;
    w_L=weg;
    Disp=disp_pna;
    
    % Compute E(3):E(5) prefactor
    alpha2=polarizability(K,weg_Solv,w_L);
    prefactor_OffRes=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t)/4/pi^2;
    prefactor_Res=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_t);
    
    %% Consolidate parameters
    % Consolidate material parameters into structure variable
    % parameters_material
    parameters_material.disp=Disp;
    parameters_material.gamma_eg=gamma_eg;
    parameters_material.gamma_vib=gamma_vib;
    parameters_material.weg=weg;
    parameters_material.wvib=wvib;
    parameters_material.wsolv=wsolv;
    
    %     % Consolidate laser parameters into structure variable parameters_laser
    parameters_laser.dt=dt;
    parameters_laser.nt=Full_2dRR_nt;
    parameters_laser.w_L=w_L;
    
    %% Compute matrix of basis states, vibrational energies, and overlap integrals
    % Basis states and vibrational energies
    [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
    % Overlap integrals
    [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);
    
    [Full_2dRR_OffRes_Cascade(:,:,iv),Full_2dRR_OffRes_Direct(:,:,iv)] = cascade_2dRR_OffResSpec(...
        wviball_TC,nquanta,ovlp_TC,parameters_laser,parameters_material);
    Full_2dRR_OffRes_Cascade(:,:,iv)=prefactor_OffRes*3e10*Full_2dRR_OffRes_Cascade(:,:,iv);
    [Full_2dRR_Res_Cascade(:,:,iv),Full_2dRR_Res_Direct(:,:,iv)] = cascade_2dRR_ResSpec(...
        wviball_TC,nquanta,ovlp_TC,parameters_laser,parameters_material);
    Full_2dRR_Res_Cascade(:,:,iv)=prefactor_Res/3e10*Full_2dRR_Res_Cascade(:,:,iv);
    
end
save('Full2dRR','Full_2dRR_OffRes_Cascade','Full_2dRR_OffRes_Direct',...
    'Full_2dRR_Res_Cascade','Full_2dRR_Res_Direct',...
    'Full_2dRR_dts','Full_2dRR_nt');
matlabpool close
'done'