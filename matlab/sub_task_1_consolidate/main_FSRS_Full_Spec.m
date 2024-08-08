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
    wvib,Disp]=PNA_parameters('dichloromethane');
% Vibrational modes to be examined in cm^-1
wvibs=wvib([1 3]);
% Unitless mode displacements
disp_pna=Disp([1 3]);
% Electronic dephasing in cm^-1
gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
% Transition dipole in D
mu_eg=translen*1e-10*q_e/Dconv;

FSRS_Full_dts=10;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;
% Compute Solvent concentration in mol m^-3
C_Solv=rho_Solv/mm_Solv;
C_Sol=2e-4;

%% Default laser paramters
% Spectral width actinic pump
LAMBDA_ap=210;
% Spectral width raman pump
LAMBDA_rp=40;
% Frequency actinic pump
w_ap=weg;
% Frequency raman pump
w_rp=w_ap;
%
FSRS_Full_nt=128;

FSRS_Full_OffRes_Cascade=complex(zeros(FSRS_Full_nt,FSRS_Full_nt,2,'double'));
FSRS_Full_OffRes_Direct=complex(zeros(FSRS_Full_nt,FSRS_Full_nt,2,'double'));
FSRS_Full_Res_Cascade=complex(zeros(FSRS_Full_nt,FSRS_Full_nt,2,'double'));
FSRS_Full_Res_Direct=complex(zeros(FSRS_Full_nt,FSRS_Full_nt,2,'double'));
matlabpool open 10
for iv=1:2
    wvib=wvibs(iv);
    dt=FSRS_Full_dts;
    Disp=disp_pna;
    w_ap=weg;
    w_rp=w_ap;
    % Compute E(3):E(5) prefactor
    alpha2=polarizability(K,weg_Solv,w_ap);
    prefactor_OffRes=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_ap-wvib)/4/pi^2;
    prefactor_Res=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_ap-wvib);
    
    %% Compute matrix of basis states, vibrational energies, and overlap integrals
    % Basis states and vibrational energies
    [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
    % Overlap integrals
    [ovlp_TC] = fcinfo_TC(base_TC,disp_pna,nmode,nquanta);
    
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
    parameters_laser.nt=FSRS_Full_nt;
    
    %% Compute third order cascade response functions
    [FSRS_Full_OffRes_Cascade(:,:,iv),FSRS_Full_OffRes_Direct(:,:,iv)] = cascade_FSRS_OffResSpec(...
        wviball_TC,2*nquanta,ovlp_TC,parameters_laser,parameters_material);
    FSRS_Full_OffRes_Cascade(:,:,iv)=prefactor_OffRes*3e10*FSRS_Full_OffRes_Cascade(:,:,iv);
    
    [FSRS_Full_Res_Cascade(:,:,iv),FSRS_Full_Res_Direct(:,:,iv)] = cascade_FSRS_ResSpec(...
        wviball_TC,2*nquanta,ovlp_TC,parameters_laser,parameters_material);
    FSRS_Full_Res_Cascade(:,:,iv)=prefactor_Res/3e10*FSRS_Full_Res_Cascade(:,:,iv);
end
save('FSRSFull','FSRS_Full_OffRes_Cascade','FSRS_Full_OffRes_Direct',...
    'FSRS_Full_Res_Cascade','FSRS_Full_Res_Direct',...
    'FSRS_Full_dts','FSRS_Full_nt');
matlabpool close
'done'