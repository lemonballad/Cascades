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
    wvibs=wvibs([1 5]);disp_pna=disp([1 5]);dts=1./wvibs/4/c;
end
% Set number of vibrational modes
nmode=1;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=3;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);

for iv=1:2
    %% Default laser paramters
    wvib=wvibs(iv);
    dt=dts(iv);
    w_L=weg;
    nt=100;
    %
    dws=000:200:2000;nw=length(dws);
    disps=(0.25:0.25:3);nd=length(disps);
   for iw=1:nw
        for id=1:nd
            w_L=dws(iw)+weg-wvib;
            disp=disps(id);
            
            % Compute E(3):E(5) prefactor
            prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_L-wvib);
                        
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);
            tic
            [rat,cas,dir] = cascadesub4_TCXi(wviball_TC,nquanta,ovlp_TC,...
                disp,gamma_eg,gamma_vib,weg,dt,nt,w_L,wvib);
            ratio(id,iw,iv)=rat*prefactor/3e10;rcas(id,iw,iv)=cas;rdir(id,iw,iv)=dir;
            [toc iv iw id]
        end
    end
end