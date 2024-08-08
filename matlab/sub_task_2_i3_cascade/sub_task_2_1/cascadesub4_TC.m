function [ratio] = cascadesub4_TC(E_vib,nquanta,ovlp,parameters_material)
% Compute third order cascade response functions
%   wviball             : vibrational energies for each basis state
%   ovlp                  : overlap integrals between basis states.
%   nquanta           : number of quanta
%   parameters      : Material parameters

%% COMPUTE BOLTZMANN POPULATIONS
kT=200; % kT in cm^-1
% Compute array of 
boltz_factor=exp(-E_vib/kT);
% Compute partition function
partition_func=sum(boltz_factor);
% Compute boltzman populations
boltz_pop=boltz_factor/partition_func;

%% Material parameters
% Electronic dephasing rate
gamma_eg=parameters_material.gamma_eg;
% Vibrational dephasing rate
gamma_vib=parameters_material.gamma_vib;
% Electronic energy gap
weg=parameters_material.weg;
% Vibrational frequencies
wvib=parameters_material.wvib;

%%
iq=nquanta;
% Speed of light in cm/fs
c=3e-5;%2.998E-5;%
% 2pi*c
r2w=2*pi*c;%0.0001885;%
% Vibrational energy gaps
[col,row]=meshgrid(E_vib,E_vib);
w=row-col;

dt=100;
wL=weg+2810;
%
kj=1:200;
tau_1=(kj-1)*dt+0;
tau_2=(kj-1)*dt+0;
nt=length(kj);

% Initialize terms of response function to 0
r1=complex(zeros(nt,nt,'double'));
r2=r1;r3=r1;r4=r1;
f1=complex(zeros(1,nt,'double'));
f2=f1;

for m=1:1%iq   
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
f1=f1+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
*1/(wL-weg-w(n,m)+1i*gamma_eg)...
.*exp(-1i*w(k,m)*r2w*tau_1-km*gamma_vib*tau_1)...
*1/(wL-weg-w(l,m)+1i*gamma_eg);
%
f2=f2+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
*1/(-wL+weg-w(m,n)+1i*gamma_eg)...
.*exp(-1i*w(m,k)*r2w*tau_1-km*gamma_vib*tau_1)...
*1/(wL-weg-w(l,k)+1i*gamma_eg);
% end
%
% DIRECT FIFTH ORDER SIGNAL
%
for ii=1:length(tau_2)
for u=1:iq
for v=1:iq
if u~=m && k~=m
r1(ii,:)=r1(ii,:)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
*1/(wL-weg-w(n,m)+1i*gamma_eg)...
.*exp(-1i*w(k,m)*r2w*tau_1-gamma_vib*tau_1)...
*1/(wL-weg-w(l,m)+1i*gamma_eg)...
.*exp(-1i*w(u,m)*r2w*tau_2(ii)-gamma_vib*tau_2(ii))...
*1/(wL-weg-w(v,m)+1i*gamma_eg);
%
r3(ii,:)=r3(ii,:)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
*1/(-wL+weg-w(m,n)+1i*gamma_eg)...
.*exp(-1i*w(m,k)*r2w*tau_1-gamma_vib*tau_1)...
*1/(-wL+weg-w(m,l)+1i*gamma_eg)...
.*exp(-1i*w(m,u)*r2w*tau_2(ii)-gamma_vib*tau_2(ii))...
*1/(wL-weg-w(v,u)+1i*gamma_eg);
end  
%
if u~=k && k~=m
r2(ii,:)=r2(ii,:)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
*1/(wL-weg-w(n,m)+1i*gamma_eg)...
.*exp(-1i*w(k,m)*r2w*tau_1-gamma_vib*tau_1)...
*1/(-wL+weg-w(k,l)+1i*gamma_eg)...
.*exp(-1i*w(k,u)*r2w*tau_2(ii)-gamma_vib*tau_2(ii))...
*1/(wL-weg-w(v,u)+1i*gamma_eg);

r4(ii,:)=r4(ii,:)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
*1/(-wL+weg-w(m,n)+1i*gamma_eg)...
.*exp(-1i*w(m,k)*r2w*tau_1-gamma_vib*tau_1)...
*1/(wL-weg-w(l,k)+1i*gamma_eg)...
.*exp(-1i*w(u,k)*r2w*tau_2(ii)-gamma_vib*tau_2(ii))...
*1/(wL-weg-w(v,k)+1i*gamma_eg);
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
direct=-0.09*(r1+r2+r3+r4);
direct=direct-mean(mean(direct));
%
% COMPUTE 2D CASCADES
%
nt1=length(kj)/2;
direct=direct(1:nt1,1:nt1);
for j=1:nt1
for jj=1:nt1
seq1(jj,j)=-7.9e-4*1i*(f1(j)+f2(j))*(f1(jj)+f2(jj));
seq2(jj,j)=-2.7e-3*1i*(conj(f1(j))+conj(f2(j)))*(f1(jj)+f2(jj));
par1(jj,j)=-8e-4*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
par2(jj,j)=-0.09*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
end
end
%
cascade=seq1+seq2+par1+par2;
cascade=cascade-mean(mean(cascade));

%
% FOURIER TRANSFORM tau_1 & tau_2
%
cascade2d=fftshift(fft2(cascade));
direct2d=fftshift(fft2(direct));


cascade2d=(cascade2d)-(cascade2d(1,100));
direct2d=(direct2d)-(direct2d(1,100));
ratio=abs(cascade2d(85,85))/abs(direct2d(85,85))*3.8954e14/3e10;

% 
% ijk=1:nt1;
dw=1/nt1;
ff=(-1/2:dw:(1/2-dw))/dt/c;
% ijk=1:nt1;
% ff(ijk)=-3.141592/dt+2*3.141592*(ijk-1)/nt1/dt;
% ff(ijk)=ff(ijk)/.0001885;
% cascade2d=abs(cascade2d)-abs(cascade2d(1,100));
% direct2d=abs(direct2d)-abs(direct2d(1,100));
% subplot(2,2,1);contour(ff,ff,abs(cascade2d).^.5,20);colorbar;title('CASCADES')
% subplot(2,2,2);contour(ff,ff,abs(direct2d).^.5,20);colorbar;title('DIRECT')
%
