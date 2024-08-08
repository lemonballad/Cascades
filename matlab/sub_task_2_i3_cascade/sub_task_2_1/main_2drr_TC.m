%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
% Plot spectra from response functions
plot_flag=true

%% Default material parameters
% Concentratin in mol L^-1
C=2E-4; % just a number
% unitless displacement
disp=7;%[0.35 0.35 1];
% electronic dephasing in cm^-1
gamma_eg=2010;
% vibrational dephasing in cm^-1
gamma_vib=1/1200;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=8.8; % this is not correct, just a number
% refractive index
n_w_t=1.39; % also just a number
% Electronic energy gap origin
weg=38000;
% Vibrational frequencies
wvib=112;
% Signal frequency
w_t=weg;
% Myoglobin parameters
if false
    [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters;
end
% PNA parameters
if false
     [gamma_vib,gamma_eg,kappa,lambda,weg,translen,n_w_t,...
    wvib,disp]=PNA_parameters('acetonitrile');
    wvib=wvib(1);disp=disp(1)
end
% Set number of vibrational modes
nmode=length(wvib);
% Set number of vibrational quanta to apportion
nquanta=7;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);
  
%% Consolidate parameters
% Consolidate material parameters into structure variable
% parameters_material
parameters_material.disp=disp;
parameters_material.gamma_eg=gamma_eg;
parameters_material.gamma_vib=gamma_vib;
parameters_material.weg=weg;
parameters_material.wvib=wvib;

%% Compute matrix of basis states, vibrational energies, and overlap integrals
% Basis states and vibrational energies
[base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
% Overlap integrals
[ovlp_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);

ratio = cascadesub4_TC2(wviball_TC,nquanta,ovlp_TC,...
    parameters_material,plot_flag);
