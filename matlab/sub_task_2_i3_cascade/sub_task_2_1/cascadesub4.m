function [ratio] = cascadesub4(d,gam_eg)
% d=7
% gam_eg=4022/2;
%
%
%
dt=100;
iq=3;
wL=38000+2810;
weg=38000;
wvib=112;
% gam_eg=1000; %cm-1
gam_vib=1/1200; %1/fs
[ovlp] = fcfac2(d);
r2w=0.0001885;
kT=200;
boltz_pop(1:iq)=exp(-((1:iq)-1)*wvib/kT) / (1/(1-exp(-wvib/kT)) ); %boltzmann populations
%
%
%
kj=1:200;
t1(kj)=(kj-1)*dt+0;
t2=(kj-1)*dt+0;
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
                f1(kj)=f1(kj)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                    *1/(wL-weg-(n-m)*wvib+1i*gam_eg)...
                    .*exp(-1i*(k-m)*wvib*r2w*t1(kj)-km*gam_vib*t1(kj))...
                    *1/(wL-weg-(l-m)*wvib+1i*gam_eg);
                %
                f2(kj)=f2(kj)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                    *1/(-wL+weg-(m-n)*wvib+1i*gam_eg)...
                    .*exp(-1i*(m-k)*wvib*r2w*t1(kj)-km*gam_vib*t1(kj))...
                    *1/(wL-weg-(l-k)*wvib+1i*gam_eg);
                % end
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
                                r1(ii,kj)=r1(ii,kj)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                    *1/(wL-weg-(n-m)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(k-m)*wvib*r2w*t1(kj)-km*gam_vib*t1(kj))...
                                    *1/(wL-weg-(l-m)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(u-m)*wvib*r2w*t2(ii)-um*gam_vib*t2(ii))...
                                    *1/(wL-weg-(v-m)*wvib+1i*gam_eg);
                                % r1(ii,kj)=r1(ii,kj)+1.*exp(-1i*(k-m)*wvib*r2w*t1(kj)-gam_vib*t1(kj)).*exp(-1i*(u-m)*wvib*r2w*t2(ii)-gam_vib*t2(ii));
                                %
                                r3(ii,kj)=r3(ii,kj)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                    *1/(-wL+weg-(m-n)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(m-k)*wvib*r2w*t1(kj)-km*gam_vib*t1(kj))...
                                    *1/(-wL+weg-(m-l)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(m-u)*wvib*r2w*t2(ii)-um*gam_vib*t2(ii))...
                                    *1/(wL-weg-(v-u)*wvib+1i*gam_eg);
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
                                r2(ii,kj)=r2(ii,kj)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                    *1/(wL-weg-(n-m)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(k-m)*wvib*r2w*t1(kj)-mk*gam_vib*t1(kj))...
                                    *1/(-wL+weg-(k-l)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(k-u)*wvib*r2w*t2(ii)-uk*gam_vib*t2(ii))...
                                    *1/(wL-weg-(v-u)*wvib+1i*gam_eg);
                                
                                r4(ii,kj)=r4(ii,kj)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                    *1/(-wL+weg-(m-n)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(m-k)*wvib*r2w*t1(kj)-mk*gam_vib*t1(kj))...
                                    *1/(wL-weg-(l-k)*wvib+1i*gam_eg)...
                                    .*exp(-1i*(u-k)*wvib*r2w*t2(ii)-uk*gam_vib*t2(ii))...
                                    *1/(wL-weg-(v-k)*wvib+1i*gam_eg);
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
for j=1:nt1
    for jj=1:nt1
        seq1(jj,j)=-7.9e-4*1i*(f2(j)*f2(jj)+f2(j)*f1(jj)+f1(j)*f2(jj)+f1(j)*f1(jj));
        seq2(jj,j)=-2.7e-3*1i*( conj(f2(j))*f2(jj)+conj(f2(j).*f1(jj))+conj(f1(j))*f2(jj)+conj(f1(j))*f1(jj));
        par1(jj,j)=-8e-4*1i*(f2(j+jj)*f2(jj)+f2(j+jj)*f1(jj)+f1(j+jj)*f2(jj)+f1(j+jj)*f1(jj));
        par2(jj,j)=-0.09*1i*(f2(j+jj)*f2(jj)+f2(j+jj)*f1(jj)+f1(j+jj)*f2(jj)+f1(j+jj)*f1(jj));
    end
end
%
cascade=seq1+seq2+par1+par2;
cascade=cascade-mean(mean(cascade));
% cascade=cascade-cascade(100,100);
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
% cascade2d=abs(cascade2d)-abs(cascade2d(1,100));
% direct2d=abs(direct2d)-abs(direct2d(1,100));
cascade2d=(cascade2d)-(cascade2d(1,100));
direct2d=(direct2d)-(direct2d(1,100));
% subplot(2,2,1);contour(ff,ff,abs(cascade2d).^.5,20);colorbar;title('CASCADES')
% subplot(2,2,2);contour(ff,ff,abs(direct2d).^.5,20);colorbar;title('DIRECT')

%
%
ratio=abs(cascade2d(85,85))/abs(direct2d(85,85))*3.8954e14/3e10;