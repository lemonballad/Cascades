clear all
% for ii=1:20
c=2.998E-5;
dt=1;
nt=1000;%+ii*100;
t=0:dt:(nt-1)*dt;
dw=1/nt;
gam=20;
wv=700;
w=(-1/2:dw:(1/2-dw))/dt/c;
wr=w(end)-w(1);
[~,iw]=min(abs(w-wv));

nf=exp((-1i*(2*pi*c)*wv-(2*pi*c)*gam)*t);
nfc=conj(nf);
fnf=fftshift(fft(nf))*2*pi*dt*c;fnf=fnf-fnf(end);
fnfc=fftshift(fft(nfc))*2*pi*dt*c;fnfc=fnfc-fnfc(end);
fnfnfc=fftshift(fft2(nfc.'*nf))*(2*pi*dt*c)^2;fnfnfc=fnfnfc-fnfnfc(end,end);
cfnf=fftshift(fft(nf.*nf))*2*pi*dt*c;cfnf=cfnf-cfnf(end);
cfnfc=fftshift(fft(nfc.*nfc))*2*pi*dt*c;cfnfc=cfnfc-cfnfc(end);
nfcfnf2=fftshift(fft2(nf.'*(nf.*nf)))*(2*pi*dt*c)^2;nfcfnf2=nfcfnf2-nfcfnf2(end,end);
nfccfnfc2=fftshift(fft2(nfc.'*(nfc.*nfc)))*(2*pi*dt*c)^2;nfccfnfc2=nfccfnfc2-nfccfnfc2(end,end);
cfnf2=fnf.'*cfnf;
cfnfc2=fnfc.'*cfnfc;
af=1i./(-wv-w+1i*gam);
afc=1i./(-wv+w+1i*gam);
afcaf=afc.'*af;
caf=conv(af,af,'same')*dt*wr/nt/2/pi;
cafc=conv(afc,afc,'same')*dt*wr/nt/2/pi;
caf2=(af.')*caf;
cafc2=afc.'*cafc;

% prop(ii)=abs(afcaf(iw,iw))/abs(fnfnfc(iw,iw));
% end
% figure;plot(prop);
% figure;plot(w,abs(fnf),'b-',w,abs(af),'r-',w,abs(fnfc),'b:',w,abs(afc),'r:')
% figure;plot(w,abs(cfnf),'b-',w,abs(caf),'r-',w,abs(cfnfc),'b:',w,abs(cafc),'r:')
figure;contourf(w,w,abs(caf2),50,'edgecolor','none');colorbar;colormap('jet');xlim([-2000 2000]);ylim([-2000 2000]);
figure;contourf(w,w,abs(nfcfnf2),50,'edgecolor','none');colorbar;colormap('jet');xlim([-2000 2000]);ylim([-2000 2000]);
figure;contourf(w,w,abs(caf2)-abs(nfcfnf2),50,'edgecolor','none');colorbar;colormap('jet');xlim([-2000 2000]);ylim([-2000 2000]);



figure;contourf(abs(f),50,'edgecolor','none');colorbar;colormap('jet');
figure;contourf(abs(f2),50,'edgecolor','none');colorbar;colormap('jet');
figure;contourf(abs(f)-abs(f2),50,'edgecolor','none');colorbar;colormap('jet');