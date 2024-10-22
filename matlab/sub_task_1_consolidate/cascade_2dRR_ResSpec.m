function [cas2,dir2] = cascade_2dRR_ResSpec(E_vib,nquanta,ovlp,...
    parameters_laser,parameters_material)
% Compute third order cascade response functions
%   wviball             : vibrational energies for each basis state
%   ovlp                  : overlap integrals between basis states.
%   nquanta           : number of quanta
%   parameters      : Material parameters
%
%   seq                  : Sequential third order cascades
%   par                  : Parallel third order cascades
%   cascades2d    : Sum of sequential and parallel third order cascades
%   direct2d          : Direct fifth order signal

%% COMPUTE BOLTZMANN POPULATIONS
kT=200; % kT in cm^-1
% Compute array of 
boltz_factor=exp(-E_vib/kT);
% Compute partition function
partition_func=sum(boltz_factor);
% Compute boltzman populations
boltz_pop=boltz_factor/partition_func;

%% Laser parameters
dt=parameters_laser.dt;
nt=parameters_laser.nt;
w_L=parameters_laser.w_L;

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
c=2.998E-5;%
% Vibrational energy gaps
[wi,wf]=meshgrid(E_vib,E_vib);
w=wf-wi;

% Initialize terms of response function to 0
dw=1/nt;
ww=(-1/2:dw:(1/2-dw))/dt/c;
nw=length(ww);
mw=max(abs(1./(-wvib-E_vib+1i*gamma_vib)));

r1=complex(zeros(nw,nw,iq,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1PA=complex(zeros(nw,iq,'double'));
f2PA=f1PA;f3PA=f1PA;f4PA=f1PA;
f1PC=complex(zeros(nw,iq,'double'));
f2PC=f1PC;f3PC=f1PC;f4PC=f1PC;
f1P=complex(zeros(nw,nw,iq,'double'));
f2P=f1P;f3P=f1P;f4P=f1P;

[W1,W2,Wf,Wi]=ndgrid(ww,ww,E_vib,E_vib);
W1=-W1-Wf+Wi;W2=-W2-Wf+Wi;
[W,Wf,Wi]=ndgrid(ww,E_vib,E_vib);
WC=-W+Wf-Wi;
W=-W-Wf+Wi;
Lp=1./(w_L-weg-w+1i*gamma_eg);
Lm=1./(-w_L+weg-w+1i*gamma_eg);
Lc1=1./(W1+1i*gamma_vib);
Lc2=1./(W2+1i*gamma_vib);
LcP=1./(W+1i*gamma_vib);
LcP2d=1./(W1+1i*gamma_vib)./(W2+1i*gamma_vib)/mw;
LcPC=1./(WC+1i*gamma_vib);
%% Loop over initial state m => assume m is in the ground state m=1, 0
% quanta.
for m=1:1
    % Loop over vibrational states n
    parfor n=1:iq
        % Loop over vibrational states k
        for k=1:iq
            % Loop over vibrational states l
            for l=1:iq
                %
                % NOTE OVERLAP INEGRALS ARE STORED WITH EXCITED STATE INDEX FIRST
                %
                if k~=m
                    f1PA(:,n)=f1PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP(:,k,m)*Lp(l,m);
                    f2PA(:,n)=f2PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP(:,m,k)*Lp(l,k);
                    f1P(:,:,n)=f1P(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP2d(:,:,k,m)*Lp(l,m);
                    f2P(:,:,n)=f2P(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP2d(:,:,m,k)*Lp(l,k);
                    f1PC(:,n)=f1PC(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcPC(:,k,m)*Lp(l,m);
                    f2PC(:,n)=f2PC(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcPC(:,m,k)*Lp(l,k);
                end
                %
                if k~=n
                    f3PA(:,n)=f3PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lp(n,l);
                    f4PA(:,n)=f4PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lp(k,l);
                    f3P(:,:,n)=f3P(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP2d(:,:,n,k)*Lp(n,l);
                    f4P(:,:,n)=f4P(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP2d(:,:,k,n)*Lp(k,l);
                    f3PC(:,n)=f3PC(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcPC(:,n,k)*Lp(n,l);
                    f4PC(:,n)=f4PC(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcPC(:,k,n)*Lp(k,l);
                end
                %
                % Compute fifth order auxillary functions
                % Loop over vibrational states u
                for u=1:iq
                    % Loop over vibrational states v
                    for v=1:iq
                        %
                        if m~=u && m~=k
                            r1(:,:,n)=r1(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                *Lp(n,m)...
                                *Lc1(:,:,k,m)...
                                *Lp(l,m)...
                                .*Lc2(:,:,u,m)...
                                *Lp(v,m);
                            
                        end
                        %
                        if k~=u && m~=k
                            r2(:,:,n)=r2(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                *Lp(n,m)...
                                *Lc1(:,:,k,m)...
                                *Lm(k,l)...
                                .*Lc2(:,:,k,u)...
                                *Lp(v,u);
                        end
                        %
                        if m~=u && m~=k
                            r3(:,:,n)=r3(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                *Lm(m,n)...
                                *Lc1(:,:,m,k)...
                                *Lm(m,l)...
                                .*Lc2(:,:,m,u)...
                                *Lp(v,u);
                        end
                        %
                        if k~=u && m~=k
                            r4(:,:,n)=r4(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                *Lm(m,n)...
                                *Lc1(:,:,m,k)...
                                *Lp(l,k)...
                                .*Lc2(:,:,u,k)...
                                *Lp(v,k);
                        end
                        %
                        %
                        if true
                            if k~=u && n~=k
                                r5(:,:,n)=r5(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(:,:,u,k)...
                                    *Lp(u,v);
                            end
                            %
                            if n~=u && n~=k
                                r6(:,:,n)=r6(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(:,:,n,u)...
                                    *Lp(n,v);
                                %
                                r7(:,:,n)=r7(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(:,:,u,n)...
                                    *Lp(u,v);
                            end
                            %
                            if k~=u && n~=k
                                r8(:,:,n)=r8(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lp(k,l)...
                                    .*Lc2(:,:,k,u)...
                                    *Lp(k,v);
                            end
                            %
                            if l~=u && m~=k
                                r9(:,:,n)=r9(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,k,m)...
                                    *Lp(l,m)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(l,v);
                                r10(:,:,n)=r10(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,k,m)...
                                    *Lm(k,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(u,v);
                                %
                                r11(:,:,n)=r11(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,m,k)...
                                    *Lp(l,k)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(l,v);
                                %
                                r12(:,:,n)=r12(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,m,k)...
                                    *Lm(m,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(u,v);
                            end
                            if l~=u && n~=k
                                %
                                r13(:,:,n)=r13(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(v,u);
                                %
                                % Auxillary function #14 for fifth order signal
                                r14(:,:,n)=r14(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(v,l);
                                %
                                % Auxillary function #15 for fifth order signal
                                r15(:,:,n)=r15(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(v,u);
                                %
                                r16(:,:,n)=r16(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lp(k,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(v,l);
                            end
                        end
                        %
                    end % End loop over vibrational states v
                end % End loop over vibrational states u
            end % End loop over vibrational state l
        end % End loop over vibrational states k
    end % End loop over vibrational states n
end % End loop over initial state
    f1PA=sum(f1PA,2);f2PA=sum(f2PA,2);f3PA=sum(f3PA,2);f4PA=sum(f4PA,2);
    f1PC=sum(f1PC,2);f2PC=sum(f2PC,2);f3PC=sum(f3PC,2);f4PC=sum(f4PC,2);
    f1P=sum(f1P,3);f2P=sum(f2P,3);f3P=sum(f3P,3);f4P=sum(f4P,3);
    r1=sum(r1,3);r2=sum(r2,3);r3=sum(r3,3);r4=sum(r4,3);
    r5=sum(r5,3);r6=sum(r6,3);r7=sum(r7,3);r8=sum(r8,3);
    r9=sum(r9,3);r10=sum(r10,3);r11=sum(r11,3);r12=sum(r12,3);
    r13=sum(r13,3);r14=sum(r14,3);r15=sum(r15,3);r16=sum(r16,3);
%
[~,iomega]=min(abs(ww-wvib));

dir2=((r1+r2+r3+r4)+(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)+(r13+r14+r15+r16));
direct=dir2(iomega,iomega);
%
% COMPUTE 2D CASCADES
%
seq1=(f1PA+f2PA+f3PA+f4PA)*(f1PA+f2PA+f3PA+f4PA).';
seq2=(f1PC+f2PC+f3PC+f4PC)*(f1PA+f2PA+f3PA+f4PA).';
fPA=(f1PA+f2PA+f3PA+f4PA);
fP=(f1P+f2P+f3P+f4P);
wr=1/nw;%2*pi*10/wvib/nw/c;%(ww(end)-ww(1))/nw/dt;% par1=par1*2*wr/nt/c;
parfor ii=1:nw
    par1(:,ii)=conv(fPA,fP(:,ii),'same')*wr;
end
par2=par1;
seq=seq1+seq2;
par=par1+par2;
cas2=seq+par;

% figure;
% subplot(2,2,1);contour(ww,ww,abs(1.8841e+13/3e10*cas2),50);colorbar;colormap('jet');
% subplot(2,2,2);contour(ww,ww,abs(dir2),50);colorbar;colormap('jet');
% subplot(2,2,3);contour(ww,ww,abs(dir2+1.8841e+13/3e10*cas2),50);colorbar;colormap('jet');
% subplot(2,2,4);contour(ww,ww,log10(abs(1.8841e+13/3e10*cas2)./abs(dir2)),50);colorbar;colormap('jet');