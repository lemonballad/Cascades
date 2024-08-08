function [cas2,direct] = cascade_2dRR_OffResSpec(E_vib,nquanta,ovlp,...
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
% Vibrational frequencies Solute
wvib=parameters_material.wvib;
% Vibrational frequencies Solvent
wSolv=parameters_material.wsolv;
% Vibrational energy levels Solvent
wsolv=(0:nquanta+1/2)*wSolv;

%% COMPUTE BOLTZMANN POPULATIONS
kT=200; % kT in cm^-1
% Compute array of 
boltz_factor=exp(-E_vib/kT);
boltz_factor_Solv=exp(-wsolv/kT);
% Compute partition function
partition_func=sum(boltz_factor);
partition_func_Solv=sum(boltz_factor_Solv);
% Compute boltzman populations
boltz_pop=boltz_factor/partition_func;
boltz_pop_Solv=boltz_factor_Solv/partition_func_Solv;

%%
iq=nquanta;
% Speed of light in cm/fs
c=2.998E-5;%

% Vibrational energy gaps Solute
[wi,wf]=meshgrid(E_vib,E_vib);
w=wf-wi;

dw=1/nt;
ww=(-1/2:dw:(1/2-dw))/dt/c;
nw=length(ww);
mw=max(abs(1./(-wvib-E_vib+1i*gamma_vib)));
mws=max(abs(1./(-wvib-wsolv+1i*gamma_vib)));

% Initialize terms of response function to 0
r1=complex(zeros(nw,nw,iq,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1PA=complex(zeros(nw,iq,'double'));f2PA=f1PA;f3PA=f1PA;f4PA=f1PA;
f1PAC=complex(zeros(nw,iq,'double'));f2PAC=f1PAC;f3PAC=f1PAC;f4PAC=f1PAC;
f1PB=complex(zeros(nw,nw,iq,'double'));f2PB=f1PB;f3PB=f1PB;f4PB=f1PB;
f1PSA=complex(zeros(nw,iq,'double'));f2PSA=f1PSA;
f1PSAC=complex(zeros(nw,iq,'double'));f2PSAC=f1PSAC;
f1PSB=complex(zeros(nw,nw,iq,'double'));f2PSB=f1PSB;

[W1,W2,Wf,Wi]=ndgrid(ww,ww,E_vib,E_vib);
W1=-W1-Wf+Wi;W2=-W2-Wf+Wi;
[WS1,WS2,WSf,WSi]=ndgrid(ww,ww,wsolv,wsolv);
WS1=-WS1-WSf+WSi;WS2=-WS2-WSf+WSi;
[W,Wf,Wi]=ndgrid(ww,E_vib,E_vib);W=-W+Wf-Wi;
[WS,WSf,WSi]=ndgrid(ww,wsolv,wsolv);WS=-WS+WSf-WSi;

Lp=1./(w_L-weg-w+1i*gamma_eg);
Lm=1./(-w_L+weg-w+1i*gamma_eg);
Lc1=1./(W1+1i*gamma_vib);
Lc2=1./(W2+1i*gamma_vib);
LcP=1./(W+1i*gamma_vib);
LcPS=1./(WS+1i*gamma_vib);
LcPC=1./(W+1i*gamma_vib);
LcPSC=1./(WS+1i*gamma_vib);
LcP2d=1./(W1+1i*gamma_vib)./(W2+1i*gamma_vib)/mw;
LcPS2d=1./(WS1+1i*gamma_vib)./(WS2+1i*gamma_vib)/mws;
%% Loop over initial state m => assume m is in the ground state m=1, 0
% quanta.
parfor m=1:iq
    % Loop over vibrational states n
    for n=1:iq
        if abs(n-m)==1
            %
            f1PSA(:,m)=f1PSA(:,m)+boltz_pop_Solv(m)*LcPS(:,n,m);
            f2PSA(:,m)=f2PSA(:,m)+boltz_pop_Solv(m)*LcPS(:,m,n);
            f1PSAC(:,m)=f1PSAC(:,m)+boltz_pop_Solv(m)*LcPSC(:,n,m);
            f2PSAC(:,m)=f2PSAC(:,m)+boltz_pop_Solv(m)*LcPSC(:,m,n);
            f1PSB(:,:,m)=f1PSB(:,:,m)+boltz_pop_Solv(m)*LcPS2d(:,:,n,m);
            f2PSB(:,:,m)=f2PSB(:,:,m)+boltz_pop_Solv(m)*LcPS2d(:,:,m,n);
            %
        end
        % Loop over vibrational states k
        for k=1:iq
            % Loop over vibrational states l
            for l=1:iq
                %
                % NOTE OVERLAP INEGRALS ARE STORED WITH EXCITED STATE INDEX FIRST
                %
                if k~=m
                    f1PA(:,m)=f1PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP(:,k,m)*Lp(l,m);
                    f2PA(:,m)=f2PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP(:,m,k)*Lp(l,k);
                    f1PAC(:,m)=f1PAC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcPC(:,k,m)*Lp(l,m);
                    f2PAC(:,m)=f2PAC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcPC(:,m,k)*Lp(l,k);
                    f1PB(:,:,m)=f1PB(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP2d(:,:,k,m)*Lp(l,m);
                    f2PB(:,:,m)=f2PB(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP2d(:,:,m,k)*Lp(l,k);
                    %
                end
                %
                if k~=n
                    f3PA(:,m)=f3PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lp(n,l);
                    f4PA(:,m)=f4PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lp(k,l);
                    f3PAC(:,m)=f3PAC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcPC(:,n,k)*Lp(n,l);
                    f4PAC(:,m)=f4PAC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcPC(:,k,n)*Lp(k,l);
                    f3PB(:,:,m)=f3PB(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP2d(:,:,n,k)*Lp(n,l);
                    f4PB(:,:,m)=f4PB(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP2d(:,:,k,n)*Lp(k,l);
                end
                %
                % Compute fifth order auxillary functions
                % Loop over vibrational states u
                for u=1:iq
                    % Loop over vibrational states v
                    for v=1:iq
                        %
                        if m~=u && m~=k
                            r1(:,:,m)=r1(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                *Lp(n,m)...
                                *Lc1(:,:,k,m)...
                                *Lp(l,m)...
                                .*Lc2(:,:,u,m)...
                                *Lp(v,m);
                            
                        end
                        %
                        if k~=u && m~=k
                            r2(:,:,m)=r2(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                *Lp(n,m)...
                                *Lc1(:,:,k,m)...
                                *Lm(k,l)...
                                .*Lc2(:,:,k,u)...
                                *Lp(v,u);
                        end
                        %
                        if m~=u && m~=k
                            r3(:,:,m)=r3(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                *Lm(m,n)...
                                *Lc1(:,:,m,k)...
                                *Lm(m,l)...
                                .*Lc2(:,:,m,u)...
                                *Lp(v,u);
                        end
                        %
                        if k~=u && m~=k
                            r4(:,:,m)=r4(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
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
                                r5(:,:,m)=r5(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(:,:,u,k)...
                                    *Lp(u,v);
                            end
                            %
                            if n~=u && n~=k
                                r6(:,:,m)=r6(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(:,:,n,u)...
                                    *Lp(n,v);
                                %
                                r7(:,:,m)=r7(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(:,:,u,n)...
                                    *Lp(u,v);
                            end
                            %
                            if k~=u && n~=k
                                r8(:,:,m)=r8(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lp(k,l)...
                                    .*Lc2(:,:,k,u)...
                                    *Lp(k,v);
                            end
                            %
                            if l~=u && m~=k
                                r9(:,:,m)=r9(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,k,m)...
                                    *Lp(l,m)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(l,v);
                                r10(:,:,m)=r10(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,k,m)...
                                    *Lm(k,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(u,v);
                                %
                                r11(:,:,m)=r11(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,m,k)...
                                    *Lp(l,k)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(l,v);
                                %
                                r12(:,:,m)=r12(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,m,k)...
                                    *Lm(m,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(u,v);
                            end
                            if l~=u && n~=k
                                %
                                r13(:,:,m)=r13(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(v,u);
                                %
                                % Auxillary function #14 for fifth order signal
                                r14(:,:,m)=r14(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *Lp(n,m)...
                                    *Lc1(:,:,n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(:,:,u,l)...
                                    *Lp(v,l);
                                %
                                % Auxillary function #15 for fifth order signal
                                r15(:,:,m)=r15(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lm(m,n)...
                                    *Lc1(:,:,k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(:,:,l,u)...
                                    *Lp(v,u);
                                %
                                r16(:,:,m)=r16(:,:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
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
f1PAC=sum(f1PAC,2);f2PAC=sum(f2PAC,2);f3PAC=sum(f3PAC,2);f4PAC=sum(f4PAC,2);
f1PB=sum(f1PB,3);f2PB=sum(f2PB,3);f3PB=sum(f3PB,3);f4PB=sum(f4PB,3);
f1PSA=sum(f1PSA,2);f2PSA=sum(f2PSA,2);
f1PSAC=sum(f1PSAC,2);f2PSAC=sum(f2PSAC,2);
f1PSB=sum(f1PSB,3);f2PSB=sum(f2PSB,3);
r1=sum(r1,3);r2=sum(r2,3);r3=sum(r3,3);r4=sum(r4,3);
r5=sum(r5,3);r6=sum(r6,3);r7=sum(r7,3);r8=sum(r8,3);
r9=sum(r9,3);r10=sum(r10,3);r11=sum(r11,3);r12=sum(r12,3);
r13=sum(r13,3);r14=sum(r14,3);r15=sum(r15,3);r16=sum(r16,3);

%
direct=((r1+r2+r3+r4)+(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)+(r13+r14+r15+r16));

% COMPUTE 2D CASCADES
seq1=-(f1PSA+f2PSA+f1PSAC+f2PSAC)...
    *(f1PA+f2PA+f3PA+f4PA).';
seq2=-(f1PA+f2PA+f3PA+f4PA+f1PAC+f2PAC+f3PAC+f4PAC)...
    *(f1PSA+f2PSA).';
fPA=(f1PA+f2PA+f3PA+f4PA);
fPSA=(f1PSA+f2PSA);
fPB=(f1PB+f2PB+f3PB+f4PB);
fPSB=(f1PSB+f2PSB);
wr=1/nw;
parfor ii=1:nw
    par1(:,ii)=2*conv(fPA,fPSB(:,ii),'same')*wr;
    par2(:,ii)=2*conv(fPSA,fPB(:,ii),'same')*wr;
end

seq=seq1+seq2;
par=par1+par2;
cas2=seq+par;

% figure;
% subplot(2,2,1);contour(ww,ww,abs(cas2),50);colorbar;colormap('jet');
% subplot(2,2,2);contour(ww,ww,abs(direct),50);colorbar;colormap('jet');
% subplot(2,2,3);contour(ww,ww,abs(6.8e-16*3e10*cas2+direct),50);colorbar;colormap('jet');
% subplot(2,2,4);contour(ww,ww,abs(6.8e-16*3e10*cas2)+abs(direct),50);colorbar;colormap('jet');

