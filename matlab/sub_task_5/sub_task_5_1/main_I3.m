%% This script calculates third order cascades
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all

%% Set parameters for calculation
nmode=3;
nquanta=7;
wvib=[1370 1100 1500];
disp=[0.35 0.35 1];

%% Compute matrix of basis states, vibrational energies, and overlap integrals
tic
[base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
[fcall_tc] = fcinfo_TC(base_TC,disp,nmode,nquanta);
toc
%fcall=abs(fcall);

%% Compute third order cascade response functions
[direct,seq,e3,tau,wrs] = response2_TC(wviball_TC,fcall_tc);

%% Plot
%
%
fac=3e13/3e10/2;%/2/pi;
%
figure;
subplot(2,2,1);plot(wrs,real(seq(100,:))*fac,wrs,real(direct(25,:))    );
% legend('cascade','direct')
subplot(2,2,2);plot(wrs,imag(seq(100,:))*fac,wrs,imag(direct(100,:))    );
% legend('cascade','direct')
subplot(2,2,3);semilogy(wrs,abs(seq(100,:))*fac,wrs,abs(direct(100,:))    );
% legend('cascade','direct')
subplot(2,2,4);plot(wrs,abs(seq(100,:))*fac,wrs,abs(direct(100,:))    );
% figure;plot(wrs,abs(seq(100,:))./abs(direct(100,:))*2e3    );
% plot(wrs,real(direct(100,:)),wrs,real(seq(100,:))*2e3  )
figure;
subplot(2,1,1);plot(wrs,real(seq(100,:))./max(abs(real(seq(100,:)))),wrs,real(direct(100,:))./max(abs(real(direct(100,:)))),wrs,real(e3)/max(abs(real(e3)))   );
legend('cascade','direct','4wm')
subplot(2,1,2);plot(wrs,imag(seq(100,:))./max(abs(imag(seq(100,:)))),wrs,imag(direct(100,:))./max(abs(imag(direct(100,:)))),wrs,imag(e3)/max(abs(imag(e3)))   );