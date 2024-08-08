%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
% Define constants
c=2.998E-5;% Speed of light cm/fs
q_e=1.60218e-19; % Fundamental charge of an electron in C
Dconv=3.33564e-30; % Conversion factor to D
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
% unitless displacement
Disp=0.35;%[0.35 0.35 1];
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
CH2Cl2_2DRR_Vib_nt=128;
CH2Cl2_2DRR_Vib_wsolbs=0:50:2000;nw=length(CH2Cl2_2DRR_Vib_wsolbs);
CH2Cl2_2DRR_Vib_Ratio=complex(zeros(nw,2,'double'));
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
CH2Cl2_2DRR_Vib_dts=10;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;

% Compute Solvent concentration in mol m^-3
C_Solv=rho_Solv/mm_Solv;
C_Sol=2e-4;

    %% Default laser paramters
    wvib=wvibs(1);
    dt=CH2Cl2_2DRR_Vib_dts;
   for iw=1:nw
       tic

       wsolv=CH2Cl2_2DRR_Vib_wsolbs(iw);
            w_L=weg;
            Disp=disp_pna;
            
            % Compute E(3):E(5) prefactor
            alpha2=polarizability(K,weg_Solv,w_L);
            prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t)/4/pi^2;
            prefactorres=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_t);

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
            parameters_laser.nt=CH2Cl2_2DRR_Vib_nt;
            parameters_laser.w_L=w_L;
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);
            
            [rat,cas,dir] = cascade_2dRR_OffRes(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            CH2Cl2_2DRR_Vib_Ratio(iw,iv)=rat*prefactor*3e10;
            CH2Cl2_2DRR_Vib_Cas(iw,iv)=cas*prefactor*3e10;CH2Cl2_2DRR_Vib_Dir(iw,iv)=dir;
            [rat,cas,dir] = cascade_2dRR_Res(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            CH2Cl2_2DRR_Vib_Ratio2(iw,iv)=rat*prefactorres/3e10;
            CH2Cl2_2DRR_Vib_Cas2(iw,iv)=cas*prefactorres/3e10;CH2Cl2_2DRR_Vib_Dir2(iw,iv)=dir;
            

            toc
    end
end
save('CH2Cl2_2DRR_Vib','CH2Cl2_2DRR_Vib_Ratio','CH2Cl2_2DRR_Vib_Cas','CH2Cl2_2DRR_Vib_Dir'...
    ,'CH2Cl2_2DRR_Vib_Cas2','CH2Cl2_2DRR_Vib_Dir2');
'done'

figure;plot(CH2Cl2_2DRR_Vib_wsolbs,abs(CH2Cl2_2DRR_Vib_Cas(:,1)+...
    CH2Cl2_2DRR_Vib_Cas2(:,1))./abs(CH2Cl2_2DRR_Vib_Dir(:,1)),'r-','Linewidth',2);hold on
plot(CH2Cl2_2DRR_Vib_wsolbs,abs(CH2Cl2_2DRR_Vib_Cas(:,2)+...
    CH2Cl2_2DRR_Vib_Cas2(:,2))./abs(CH2Cl2_2DRR_Vib_Dir(:,2)),'b-','Linewidth',2);
xlabel('\omega_{VIB}/2\pic (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
legend('Dichloromethane','Acetonitrile');
hold off
