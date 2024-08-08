%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
c=2.998E-5;%
q_e=1.60218e-19;
Dconv=3.33564e-30;
kT=200; % kT in cm^-1

%% Default material parameters
% Solvent Raman differential cross section prefactor in m^2
K=7.62e-31; % Methanol
% Solvent density in g/m^3
rho_Solv=792000; % Methanol
% Solvent molar mass in g/mol
mm_Solv=32.04; % Methanol
% Solvent vibrational mode in cm^-1
wsolv=1035; % Methanol
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% refractive index Solute
n_w_t_Sol=1.39; % also just a number
% refractive index Solvent
n_w_t_Solv=1.34; % Methanol
% Electronic energy gap origin Solute in cm^-1
weg=22720;
% Electronic energy gap origin Solvent in cm^-1
weg_Solv=153100; % Methanol
% Signal frequency
w_t=weg;
FSRS_CH2Cl2_dws=-3000:50:3000;
nw=length(FSRS_CH2Cl2_dws);
FSRS_CH2Cl2_Ratio=complex(zeros(nw,2,'double'));
% matlabpool open 10
for iv=1:2
    if iv==1
        solvent='dichloromethane';
    else
        solvent='acetonitrile';
    end
% PNA parameters
% [~,~,kappa,lambda,weg,translen,n_w_t,...
%     wvib,Disp]=PNA_parameters('methanol');
[~,~,kappa,lambda,weg,translen,n_w_t_Sol,...
    wvib,disp,K,rho_Solv,mm_Solv,wsolv,n_w_t_Solv,weg_Solv]...
    =DEFENSE_PNA_parameters(solvent);
% Vibrational modes to be examined in cm^-1
wvibs=wvib([1 3]);
% Unitless mode displacements
disp_pna=disp([1 3]);
% Electronic dephasing in cm^-1
gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
% Transition dipole in D
mu_eg=translen*1e-10*q_e/Dconv;
% Time step for specific vibrational mode.
FSRS_CH2Cl2_dts=10;

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
%
FSRS_CH2Cl2_nt=128;

    wvib=wvibs(1);
    dt=FSRS_CH2Cl2_dts;(iv);
    for iw=1:length(FSRS_CH2Cl2_dws)
        w_ap=FSRS_CH2Cl2_dws(iw)+weg;
        w_rp=w_ap;
        tic
            Disp=disp_pna;
            % Compute E(3):E(5) prefactor
            alpha2=polarizability(K,weg_Solv,w_ap-wvib);
            prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_ap-wvib)/4/pi^2;
            prefactorres=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_ap-wvib);
            
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
            parameters_laser.nt=FSRS_CH2Cl2_nt;
            
            %% Compute third order cascade response functions
            [rat,cas,dir] = cascade_FSRS_OffRes(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            FSRS_CH2Cl2_Ratio(iw,iv)=prefactor*3e10*rat;
            rcas(iw,iv)=prefactor*3e10*cas;rdir(iw,iv)=dir;
            [rat,cas,dir] = cascade_FSRS_Res(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            rcas2(iw,iv)=prefactorres/3e10*cas;rdir2(iw,iv)=dir;
            toc
    end
end
save('FSRSCH2Cl2','FSRS_CH2Cl2_Ratio','FSRS_CH2Cl2_dws','rcas','rdir','rcas2','rdir2');
'done'

load FSRSCH2Cl2
figure;plot(FSRS_CH2Cl2_dws,abs(rcas(:,1)+...
    rcas2(:,1))./abs(rdir(:,1)),'r-','Linewidth',2);hold on
plot(FSRS_CH2Cl2_dws,abs(rcas(:,2)+...
    rcas2(:,2))./abs(rdir(:,2)),'b-','Linewidth',2);
xlabel('(\omega_{LASER}-\omega_{eg})/2\pic (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
legend('Dichloromethane','Acetonitrile');
hold off
