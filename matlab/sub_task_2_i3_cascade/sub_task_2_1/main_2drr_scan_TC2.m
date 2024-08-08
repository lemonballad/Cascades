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
gamma_eg=2010;
% vibrational dephasing in cm^-1
gamma_vib=10;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=6.2; % this is not correct, just a number
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
    wvib=wvib(1);disp=disp(1);
end
% PNA parameters
if true
     [~,gamma_eg,kappa,lambda,weg,translen,n_w_t,...
    wvib,disp]=PNA_parameters('acetonitrile');
    wvib=wvib(5);disp=disp(5);
end
% Set number of vibrational modes
nmode=length(wvib);
% Set number of vibrational quanta to apportion
nquanta=3;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);

%% Default laser paramters
dt=8.25;
w_L=weg;
nt=500;
dw=1/nt;
w=(-1/2:dw:(1/2-dw))/dt/c;
[~,iomega]=min(abs(w-wvib));

%
dws=0:20:0;nw=length(dws);
disps=(1:0.1:1);nd=length(disps);
for iw=1:nw
    for id=1:nd
        w_L=dws(iw)+weg;
        disp=disps(id);
 
        % Compute E(3):E(5) prefactor
        prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_L-wvib);
        
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
         [seq,par,cascade2d,direct2d,cascade,direct] = cascadesub4_TC2(wviball_TC,nquanta,ovlp_TC,...
            parameters_laser,parameters_material);
        ratio(id,iw)=abs(cascade2d(iomega,iomega))/abs(direct2d(iomega,iomega))*prefactor/3e10;%3.8954e14

        toc
    end
end
% figure;semilogy(disps,ratio(:,(nw+1)/2));
% figure;contour(dws,disps,ratio);

if false%plot_flag
    seq=seq-seq(1,end);
    par=par-par(1,end);
    cascade2d=(cascade2d)-(cascade2d(1,100));
    direct2d=(direct2d)-(direct2d(1,100));
    figure;
    subplot(2,3,1);contour(w,w,abs(seq).^.5,25);colorbar;title('SEQUENTIAL CASCADES')
    subplot(2,3,2);contour(w,w,abs(par).^.5,25);colorbar;title('PARALLEL CASCADES')
    subplot(2,3,3);contour(w,w,abs(cascade2d).^.5,25);colorbar;title('CASCADES')
    subplot(2,3,4);contour(w,w,abs(direct2d),25);colorbar;title('DIRECT')
    subplot(2,3,6);contour(w,w,abs(cascade2d+direct2d).^.5,25);colorbar;title('DIRECT+CASCADES');
    xlim([-1700 1700]);ylim([-1700 1700]);
end


