% clear;
% close all;
%
%
%
nmode=2;
wvib(1)=861;1370;
disp(1)=1.5;0.35;
wvib(2)=1328;
disp(2)=1.4;
wvib(3)=1500;
disp(3)=1;
%
%
%
[base,wviball]=basis(nmode,wvib);
[fcall] = fcinfo(base,disp,nmode);
%fcall=abs(fcall);
[direct,seq,e3,tau,wrs] = response2(wviball,fcall);
%
%
fac=3e13/3e10/2;%/2/pi;
%
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