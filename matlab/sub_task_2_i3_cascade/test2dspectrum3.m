clear;
d=7;gam_eg=2010;
%
%
%
% d=1;
dt=100*.5;
iq=3;
wL=38000+2810;
weg=38000;
wvib=112;
% gam_eg=1000; %cm-1
gam_vib=1/1200; %1/fs
[ovlp] = fcfac2(d);
r2w=0.0001885;
kT=200;
p(1:iq)=exp(-((1:iq)-1)*wvib/kT) / (1/(1-exp(-wvib/kT)) ); %boltzmann populations 
%
%
%
kj=1:200;
t1(kj)=(kj-1)*dt+0;
t2(kj)=(kj-1)*dt+0;
for ii=1:length(t2)
r1(ii,kj)=0;
r2(ii,kj)=0;
r3(ii,kj)=0;
r4(ii,kj)=0;
end
f1(kj)=0;
f2(kj)=0;
for m=1:iq   
for n=1:iq 
for k=1:iq 
for l=1:iq
km=1;
if k==m
km=0;
end
%
% NOTE OVERLAP INEGRALS ARE STORED WITH EXCITED STATE INDEX FIRST
%
f1(kj)=f1(kj)+p(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
*1/(wL-weg-(n-m)*wvib+1i*gam_eg)...
*1/(wL-weg-(l-m)*wvib+1i*gam_eg)...
.*exp(-1i*(k-m)*wvib*r2w*t1(kj)-km*gam_vib*t1(kj));
%
f2(kj)=f2(kj)+p(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
*1/(-wL+weg-(m-n)*wvib+1i*gam_eg)...
*1/(wL-weg-(l-k)*wvib+1i*gam_eg)...
.*exp(-1i*(m-k)*wvib*r2w*t1(kj)-km*gam_vib*t1(kj));
%  end
%
% DIRECT FIFTH ORDER SIGNAL
%
for ii=1:length(t2)
for u=1:iq
for v=1:iq
km=1;um=1;
if u~=m && k~=m
% if u==m
%     um=0;
% end
% if k==m
%     km=0;
% end
r1(ii,kj)=r1(ii,kj)+p(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
*1/(wL-weg-(n-m)*wvib+1i*gam_eg)...
*1/(wL-weg-(l-m)*wvib+1i*gam_eg)...
*1/(wL-weg-(v-m)*wvib+1i*gam_eg)...
.*exp(-1i*(k-m)*wvib*r2w*t1(kj)-gam_vib*t1(kj))...
.*exp(-1i*(u-m)*wvib*r2w*t2(ii)-gam_vib*t2(ii));
% r1(ii,kj)=r1(ii,kj)+1.*exp(-1i*(k-m)*wvib*r2w*t1(kj)-gam_vib*t1(kj)).*exp(-1i*(u-m)*wvib*r2w*t2(ii)-gam_vib*t2(ii));
%
r3(ii,kj)=r3(ii,kj)+p(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
*1/(-wL+weg-(m-n)*wvib+1i*gam_eg)...
*1/(-wL+weg-(m-l)*wvib+1i*gam_eg)...
*1/(wL-weg-(v-u)*wvib+1i*gam_eg)...
.*exp(-1i*(m-k)*wvib*r2w*t1(kj)-gam_vib*t1(kj))...
.*exp(-1i*(m-u)*wvib*r2w*t2(ii)-gam_vib*t2(ii));
end  
%
mk=1;
uk=1;
if u~=k && k~=m
%     uk=0;
% end
% if k==m
%     mk=0;
% end
r2(ii,kj)=r2(ii,kj)+p(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
*1/(wL-weg-(n-m)*wvib+1i*gam_eg)...
*1/(-wL+weg-(k-l)*wvib+1i*gam_eg)...
*1/(wL-weg-(v-u)*wvib+1i*gam_eg)...
.*exp(-1i*(k-m)*wvib*r2w*t1(kj)-gam_vib*t1(kj))...
.*exp(-1i*(k-u)*wvib*r2w*t2(ii)-gam_vib*t2(ii));

r4(ii,kj)=r4(ii,kj)+p(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
*1/(-wL+weg-(m-n)*wvib+1i*gam_eg)...
*1/(wL-weg-(l-k)*wvib+1i*gam_eg)...
*1/(wL-weg-(v-k)*wvib+1i*gam_eg)...
.*exp(-1i*(m-k)*wvib*r2w*t1(kj)-gam_vib*t1(kj))...
.*exp(-1i*(u-k)*wvib*r2w*t2(ii)-gam_vib*t2(ii));
end
%
end
end
end
%
%
end
% end
end
end
end
%
direct=-0.09*(r1+r2+r3+r4);fd=sum(sum(abs(direct)));
direct=direct-mean(mean(direct));
%
% COMPUTE 2D CASCADES
%
nt1=length(kj)/2;
for j=1:nt1
for jj=1:nt1
seq1(jj,j)=-7.9e-4*(f2(j)*f2(jj)+f2(j)*f1(jj)+f1(j)*f2(jj)+f1(j)*f1(jj));
seq2(jj,j)=-2.7e-3*( conj(f2(j))*f2(jj)+conj(f2(j).*f1(jj))+conj(f1(j))*f2(jj)+conj(f1(j))*f1(jj));
par1(jj,j)=-8e-4*(f2(j+jj)*f2(jj)+f2(j+jj)*f1(jj)+f1(j+jj)*f2(jj)+f1(j+jj)*f1(jj));
par2(jj,j)=-0.09*(f2(j+jj)*f2(jj)+f2(j+jj)*f1(jj)+f1(j+jj)*f2(jj)+f1(j+jj)*f1(jj));
end
end
%
cascade=seq1+seq2+par1+par2;
fc=sum(sum(abs(cascade)));
% cascade=seq1/7.9e-4+0*seq2/2.7e-3;
% cascade=cascade-mean(mean(cascade));
%
% SUBTRACT NON-OSCILLATORY PART OF CASCADE
%
[ii,jj]=size(cascade);
kj=1:jj;
num=1;
kj=1:ii;
for k=1:ii
cascade2(k,kj)=cascade(k,kj)-mean(cascade(k,num:ii));
end
for k=1:ii
cascade(kj,k)=cascade2(kj,k)-mean(cascade2(num:ii,k));
end
cascade=cascade-cascade(100,100);
%
% FOURIER TRANSFORM t2
%
for j=1:nt1
x0(1:nt1)=cascade(j,1:nt1);
z0(1:nt1)=direct(j,1:nt1);
x1=fft(x0);
z1=fft(z0);
for ij=1:1:length(x0)/2
x2(nt1/2+ij)=x1(ij);
x2(ij)=x1(nt1/2+ij);
z2(nt1/2+ij)=z1(ij);
z2(ij)=z1(nt1/2+ij);
end
x3(j,1:nt1)=x2(1:nt1);
z3(j,1:nt1)=z2(1:nt1);
end
%
% FOURIER TRANSFORM t1
%
x0=0;x1=0;x2=0;z0=0;z1=0;z2=0;
% x3=x3-mean(mean(x3));
for j=1:nt1
x0(1:nt1)=x3(1:nt1,j);
x1=fft(x0);
z0(1:nt1)=z3(1:nt1,j);
z1=fft(z0);
for ij=1:1:length(x0)/2
x2(nt1/2+ij)=x1(ij);
x2(ij)=x1(nt1/2+ij);
z2(nt1/2+ij)=z1(ij);
z2(ij)=z1(nt1/2+ij);
end
cascade2d(1:nt1,j)=x2(1:nt1);
direct2d(1:nt1,j)=z2(1:nt1);
end
ijk=1:nt1;
ff(ijk)=-3.141592/dt+2*3.141592*(ijk-1)/nt1/dt;
ff(ijk)=ff(ijk)/.0001885;
% 
% cascade2d=(cascade2d)-(cascade2d(1,100));
direct2d=(direct2d)-(direct2d(1,100));
subplot(2,2,1);contour(ff,ff,abs(cascade2d).^1,20);colorbar;title('CASCADES')
subplot(2,2,2);contour(ff,ff,abs(direct2d).^1,20);colorbar;title('DIRECT')
% figure;
% subplot(2,2,1);contour(ff,ff,real(cascade2d),20);colorbar;title('CASCADES')
% subplot(2,2,2);contour(ff,ff,real(direct2d),20);colorbar;title('DIRECT')
% subplot(2,2,3);contour(ff,ff,imag(cascade2d),20);colorbar;title('CASCADES')
% subplot(2,2,4);contour(ff,ff,imag(direct2d),20);colorbar;title('DIRECT')
%
%
ratio=abs(cascade2d(85,85))/abs(direct2d(85,85))*3.8954e14/3e10
% % ratio2=abs(cascade2d(50,85))/abs(direct2d(50,85))*3.8865e14/3e10
% 
dlmwrite('.\absdirect2d.dat',abs(direct2d));
dlmwrite('.\redirect2d.dat',real(direct2d));
dlmwrite('.\imdirect2d.dat',imag(direct2d));
dlmwrite('.\abscascade2d.dat',abs(cascade2d));
dlmwrite('.\recascade2d.dat',real(cascade2d));
dlmwrite('.\imcascade2d.dat',imag(cascade2d));
dlmwrite('.\ff.dat',ff);
% 
% [ii,jj]=size(cascade);
% kj=1:jj;
% num=1;
% kj=1:ii;
% for k=1:ii
% cascade2(k,kj)=cascade(k,kj)-mean(cascade(k,num:ii));
% end
% for k=1:ii
% cascade(kj,k)=cascade2(kj,k)-mean(cascade2(num:ii,k));
% end
% figure
% contourf(cascade')