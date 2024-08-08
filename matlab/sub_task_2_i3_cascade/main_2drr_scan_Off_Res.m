%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
% Plot spectra from response functions
plot_flag=true;
% Speed of light cm/fs
c=2.998E-5;

%% Default material parameters
% Solvent Raman differential cross section prefactor in m^2
K=7.62e-31; % Methanol
% Solvent density in g/m^3
rho_Solv=792000; % Methanol
% Solvent molar mass in g/mol
mm_Solv=32.04; % Methanol
% Solvent vibrational mode in cm^-1
wsolv=861; % Methanol
% unitless displacement
disp=0.35;%[0.35 0.35 1];
% electronic dephasing in cm^-1
gamma_eg=3000;%2010;
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=6.2; % this is not correct, just a number
% refractive index Solute
n_w_t_Sol=1.39; % also just a number
% refractive index Solvent
n_w_t_Solv=1.34; % Methanol
% Electronic energy gap origin Solute in cm^-1
weg=22720;
% Electronic energy gap origin Solvent in cm^-1
weg_Solv=153100; % Methanol
% Vibrational frequencies
wvib=1600;
wvibs=[861 1599];
% Signal frequency
w_t=weg;
% Myoglobin parameters
if false
    [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters;
    wvib=wvib(1);disp=disp(1);
end
% PNA parameters
if true
     [~,gamma_eg,kappa,lambda,weg,translen,n_w_t,...
    wvibs,disp]=PNA_parameters('acetonitrile');
    wvibs=wvibs([1 5]);disp_pna=disp([1 5]);dts=1./wvibs/4/c;
end
% Set number of vibrational modes
nmode=1;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;

w_L=weg;
% Compute Solvent polarizability squared in A^4 s^8 kg^-2
alpha2=polarizability(K,weg_Solv,w_L);
% Compute Solvent concentration in mol m^-3
C_Solv=rho_Solv/mm_Solv;
C_Sol=2e-4;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t);
prefactorres=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_t);
for iv=1:2
    %% Default laser paramters
    wvib=wvibs(iv);
    dt=dts(iv);
    w_L=weg;
    nt=100;
    %
%     dws=-3000:50:3000;nw=length(dws);
%     disps=[0.05 (0.1:0.1:5)];nd=length(disps);
    dws=000:100:2000;nw=length(dws);
    disps=(0.1:0.5:3.1);nd=length(disps);
   for iw=1:nw
        for id=1:nd
            w_L=dws(iw)+weg-wvib/2;
            disp=disps(id);
            
            % Compute E(3):E(5) prefactor
            alpha2=polarizability(K,weg_Solv,w_L);
            prefactor(id,iw,iv)=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t)/4/pi^2;
            prefactorres(id,iw,iv)=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_t);

            %% Consolidate parameters
            % Consolidate material parameters into structure variable
            % parameters_material
            parameters_material.disp=disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
            parameters_material.wsolv=wsolv;
            
            %     % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.dt=dt;
            parameters_laser.nt=nt;
            parameters_laser.w_L=w_L;
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);
%             tic
%             [rat,cas,dir] = cascadesub4_TCKAPPA(wviball_TC,nquanta,ovlp_TC,...
%                 parameters_laser,parameters_material,wvib);
%             ratioK(id,iw,iv)=rat*prefactor(id,iw,iv)/3e10;rcasK(id,iw,iv)=cas;rdirK(id,iw,iv)=dir;
%             [toc]% rat cas dir]
            tic
%             [rat2,cas2,dir2] = cascade_TC_Res(wviball_TC,nquanta,ovlp_TC,...
%                 parameters_laser,parameters_material,wvib);
            [rat,cas,dir] = cascade_TC_Off_Res(wviball_TC,nquanta,ovlp_TC,...
                parameters_laser,parameters_material,wvib);
            [rat2,cas2,dir2] = cascadesub4_TCGAMMACopyCopy(wviball_TC,nquanta,ovlp_TC,...
                parameters_laser,parameters_material,wvib);
%             [rat,cas,dir] = cascade_TC_Off_Res(wviball_TC,nquanta,ovlp_TC,...
%                 parameters_laser,parameters_material,wvib);
%             ratioG(id,iw,iv)=rat*prefactorres(id,iw,iv)/3e10;rcasG(id,iw,iv)=cas;rdirG(id,iw,iv)=dir;
            ratioG2(id,iw,iv)=rat2*prefactor(id,iw,iv)*3e10;rcasG2(id,iw,iv)=cas2;rdirG2(id,iw,iv)=dir2;
            ratioG(id,iw,iv)=rat*prefactor(id,iw,iv)*3e10;rcasG(id,iw,iv)=cas;rdirG(id,iw,iv)=dir;
            [toc]% iv iw id]
        end
    end
end