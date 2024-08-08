%% This script calculates third order cascades
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all

%% Default material parameters
% Concentratin in mol L^-1
C=2E-4; % just a number
% unitless displacement
disp=0.35;%[0.35 0.35 1];
% electronic dephasing in cm^-1
gamma_eg=1000;
% vibrational dephasing in cm^-1
gamma_vib=10;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=8.8; % this is not correct, just a number
% refractive index
n_w_t=1.39; % also just a number
% Electronic energy gap origin
weg=20000;
% Vibrational frequencies
wvib=1370;%[1370 1100 1500];
% Signal frequency
w_t=weg;
% Myoglobin parameters
if false
    [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters;
end
% PNA parameters
if false
     [inhomog,gamma_0,kappa,lambda,weg,translen,refrac_index,...
    wvib,disp]=PNA_parameters('methanol');
end
% Set number of vibrational modes
nmode=length(wvib);
% Set number of vibrational quanta to apportion
nquanta=4;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);

%% Default laser paramters
% Spectral width actinic pump
LAMBDA_ap=210;
% Spectral width raman pump
LAMBDA_rp=40;
% Frequency actinic pump
w_ap=weg;
% Frequency raman pump
w_rp=weg;

%% Consolidate parameters
% Consolidate material parameters into structure variable
% parameters_material
parameters_material.disp=disp;
parameters_material.gamma_eg=gamma_eg;
parameters_material.gamma_vib=gamma_vib;
parameters_material.weg=weg;
parameters_material.wvib=wvib;

% Consolidate laser parameters into structure variable parameters_laser
parameters_laser.LAMBDA_ap=LAMBDA_ap;
parameters_laser.LAMBDA_rp=LAMBDA_rp;
parameters_laser.w_ap=w_ap;
parameters_laser.w_rp=w_rp;

%% Compute matrix of basis states, vibrational energies, and overlap integrals
tic
[base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
[fcall_tc] = fcinfo_TC(base_TC,disp,nmode,nquanta);
toc

%% Compute third order cascade response functions
% [direct,seq,e3,tau,wrs] = response2_TC(wviball_TC,fcall_tc,parameters);
[direct,cascade_3,e3,tau,wrs] = FSRS_TC(wviball_TC,fcall_tc,...
    nquanta,parameters_material,parameters_laser);

%% Plot

% 3e13 approx prefactor, 3e10 c (cm/s), ?2?
% fac=3e13/3e10/2;%/2/pi;

% adjust prefactor from s^-1 => cm^-1. This will make ratio unitless
fac=prefactor/3e10/2;

if false
figure;
subplot(2,2,1);plot(wrs,real(cascade_3(100,:))*fac,'b-',wrs,real(direct(25,:)),'r-'...
    ,'linewidth',2);
axis square;xlabel('\omega (cm^{-1})');ylabel('');xlim([0 max(wrs)]);
set(gca,'linewidth',2,'fontsize',16);title('Real Component');
subplot(2,2,2);plot(wrs,imag(cascade_3(100,:))*fac,'b-',wrs,imag(direct(100,:)),'r-'...
    ,'linewidth',2);% legend('cascade','direct')
axis square;xlabel('\omega (cm^{-1})');ylabel('');xlim([0 max(wrs)]);
set(gca,'linewidth',2,'fontsize',16);title('Imaginary Component');
subplot(2,2,3);semilogy(wrs,abs(cascade_3(100,:))*fac,'b-',wrs,abs(direct(100,:)),'r-'...
    ,'linewidth',2);% legend('cascade','direct')
axis square;xlabel('Raman Shift (cm^{-1})');xlim([0 max(wrs)]);
ylabel('Signal Field Magnitude (a.u.)');%ylim([]);
set(gca,'linewidth',2,'fontsize',16);title('log of Magnitude');
subplot(2,2,4);plot(wrs,abs(cascade_3(100,:))*fac,'b-',wrs,abs(direct(100,:)),'r-'...
    ,'linewidth',2);% legend('cascade','direct')
axis square;xlabel('\omega (cm^{-1})');ylabel('');xlim([0 max(wrs)]);
set(gca,'linewidth',2,'fontsize',16);title('Magnitude');

% figure;plot(wrs,abs(seq(100,:))./abs(direct(100,:))*2e3    );
% plot(wrs,real(direct(100,:)),wrs,real(seq(100,:))*2e3  )
figure;
subplot(1,2,1);plot(wrs,real(cascade_3(100,:))./max(abs(real(cascade_3(100,:)))),'b-',...
    wrs,real(direct(100,:))./max(abs(real(direct(100,:)))),'r-',...
    wrs,real(e3(100,:))/max(abs(real(e3(100,:)))),'g-',...
    'linewidth',2);
axis square;xlabel('\omega (cm^{-1})');ylabel('');xlim([0 max(wrs)]);
set(gca,'linewidth',2,'fontsize',16);title('Scaled Real Component');
legend('cascade','direct','4wm')
subplot(1,2,2);plot(wrs,imag(cascade_3(100,:))./max(abs(imag(cascade_3(100,:)))),'b-',...
    wrs,imag(direct(100,:))./max(abs(imag(direct(100,:)))),'r-',...
    wrs,imag(e3(100,:))/max(abs(imag(e3(100,:)))),'g-',...
    'linewidth',2);
axis square;xlabel('\omega (cm^{-1})');ylabel('');xlim([0 max(wrs)]);
set(gca,'linewidth',2,'fontsize',16);title('Scaled Imaginary Component');
end


if false
    figure;semilogy(wrs,abs(direct(100,:)),'b-',wrs,abs(cascade_3(100,:))*fac,'g-','linewidth',2    );
    ylim([1e-18 1e-13]);ylabel('Signal Field Magnitude (a.u)');
    xlim([0 3000]);xlabel('Raman Shift (cm^{-1})');
    legend('|E^{(5)}_{direct}|','|E_{cas}|');
    set(gca,'linewidth',2,'fontsize',16);axis square;caxis([1e-18 1e-13]);
end

if false   
    figure;contourf(tau,wrs,abs(direct)',50,'edgecolor','none');
    axis square;ylabel('\omega (cm^{-1})');xlabel('\tau_1');ylim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Magnitude';'of Direct Fifth Order'});
        
    figure;contourf(tau,wrs,abs(cascade_3)',50,'edgecolor','none');
    axis square;ylabel('\omega (cm^{-1})');xlabel('\tau_1');ylim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Magnitude';'of Third Order Cascades'});
        
    figure;contourf(tau,wrs,abs(cascade_3)./abs(direct),50,'edgecolor','none');
    axis square;ylabel('\omega (cm^{-1})');xlabel('\tau_1');ylim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Magnitude of Ratio';'E^{(3)}_{cascade}/E^{(5)}_{direct}'});
    
end