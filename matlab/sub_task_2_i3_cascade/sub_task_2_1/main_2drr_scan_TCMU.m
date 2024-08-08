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
% Concentratin in mol L^-1
C=2E-4; % just a number
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
% refractive index
n_w_t=1.39; % also just a number
% Electronic energy gap origin
weg=22720;
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
    wvibs=wvibs([1 5]);disp_pna=disp([1 5]);dts=1./wvibs/4/c;dts(dts==max(dts))=max(dts)/2;
end
% Set number of vibrational modes
nmode=1;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=4;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);

for iv=1:2
    %% Default laser paramters
    wvib=wvibs(iv);
    dt=7.45;%dts(iv);
    w_L=weg;
    nt=200;
    %
            w_L=+weg;
            disp=disp_pna(iv);
            
            % Compute E(3):E(5) prefactor
            prefactor(iv)=prefactor_3_5(l,C,mu_eg,n_w_t,w_L+wvib);
            
            %% Consolidate parameters
            % Consolidate material parameters into structure variable
            % parameters_material
            parameters_material.disp=disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
            
            %     % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.dt=dt;
            parameters_laser.nt=nt;
            parameters_laser.w_L=w_L;
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);
            tic
            [rat(iv),cascade2d(:,:,iv),direct2d(:,:,iv)] = cascadesub4_TCMU(wviball_TC,nquanta,ovlp_TC,...
                parameters_laser,parameters_material,wvib);
            cascade2d(:,:,iv)=prefactor(iv)*cascade2d(:,:,iv)/3e10;
%             ratio(id,iw,iv)=rat*prefactor(id,iw,iv)/3e10;rcas(id,iw,iv)=cas;rdir(id,iw,iv)=dir;
            [toc iv ]
    
end

