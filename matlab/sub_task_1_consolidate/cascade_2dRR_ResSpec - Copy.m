function [cas,dir] = cascade_TC_Res(E_vib,nquanta,ovlp,...
    parameters_laser,parameters_material,iomega)
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
% 2pi*c
r2w=2*pi*c;%0.0001885;%
% Vibrational energy gaps
[col,row]=meshgrid(E_vib,E_vib);
w=row-col;

damp=gamma_vib*c;

%
kjD=1:nt;
tauD=(kjD-1)*dt;
ntC=nt*2;
kjC=1:ntC;
tauC=(kjC-1)*dt;

% Initialize terms of response function to 0
dw=1/nt;
ww=(-1/2:dw:(1/2-dw))/dt/c;
nw=length(ww);

r1=complex(zeros(nw,nw,iq,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1PA=complex(zeros(nw,iq,'double'));
f2PA=f1PA;f3PA=f1PA;f4PA=f1PA;
f1PB=complex(zeros(nw,iq,'double'));
f2PB=f1PB;f3PB=f1PB;f4PB=f1PB;
f1PC=complex(zeros(nw,iq,'double'));
f2PC=f1PC;f3PC=f1PC;f4PC=f1PC;

[W1,W2,Wf,Wi]=ndgrid(ww,ww,E_vib,E_vib);
W1=-W1-Wf+Wi;W2=-W2-Wf+Wi;
[W,Wf,Wi]=ndgrid(ww,E_vib,E_vib);
WC=-W+Wf-Wi;
W=-W-Wf+Wi;
Lp=1i./(w_L-weg-w+1i*gamma_eg);
Lm=1i./(-w_L+weg-w+1i*gamma_eg);
Lc1=1i./(W1+1i*gamma_vib);
Lc2=1i./(W2+1i*gamma_vib);
LcP=1i./(W+1i*gamma_vib);
LcPC=1i./(WC+1i*gamma_vib);
%% Loop over initial state m => assume m is in the ground state m=1, 0
% quanta.
parfor m=1:iq
    % Loop over vibrational states n
    for n=1:iq
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
                    f1PB(:,m)=f1PB(:,m)+Lp(n,m)*LcP(:,k,m)*Lp(l,m);
                    f2PB(:,m)=f2PB(:,m)+Lm(m,n)*LcP(:,m,k)*Lp(l,k);
                    f1PC(:,m)=f1PC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcPC(:,k,m)*Lp(l,m);
                    f2PC(:,m)=f2PC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcPC(:,m,k)*Lp(l,k);
                end
                %
                if k~=n
                    f3PA(:,m)=f3PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lp(n,l);
                    f4PA(:,m)=f4PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lp(k,l);
                    f3PB(:,m)=f3PB(:,m)+Lp(n,m)*LcP(:,n,k)*Lp(n,l);
                    f4PB(:,m)=f4PB(:,m)+Lm(m,n)*LcP(:,k,n)*Lp(k,l);
                    f3PC(:,m)=f3PC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcPC(:,n,k)*Lp(n,l);
                    f4PC(:,m)=f4PC(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
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
    f1PB=sum(f1PB,2);f2PB=sum(f2PB,2);f3PB=sum(f3PB,2);f4PB=sum(f4PB,2);
    f1PC=sum(f1PC,2);f2PC=sum(f2PC,2);f3PC=sum(f3PC,2);f4PC=sum(f4PC,2);
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
wr=(ww(end)-ww(1));
seq1=(f1PA+f2PA+f3PA+f4PA)*(f1PA+f2PA+f3PA+f4PA).';
seq2=(f1PC+f2PC+f3PC+f4PC)*(f1PA+f2PA+f3PA+f4PA).';
% par1=conv((f1PA+f2PA+f3PA+f4PA),(f1PA+f2PA+f3PA+f4PA),'same')*(f1PB+f2PB+f3PB+f4PB).'*2*wr/nt/c;
par1=conv((f1PA),(f1PA),'same')*(f1PB).'*2*wr/nt/c;
par2=par1;
seq=seq1+seq2;
par=par1+par2;%par=par*2e6;
cas2=0*seq*pi/dt+par;
cascade=cas2(iomega,iomega);

cas=cas2;%abs(cascade);
dir=dir2;%abs(direct);
% ratio=cas/dir;

figure;
subplot(2,2,1);contour(ww,ww,abs(cas2),50);colorbar;colormap('jet');
subplot(2,2,2);contour(ww,ww,abs(dir2),50);colorbar;colormap('jet');
subplot(2,2,3);contour(ww,ww,abs(dir2+cas2),50);colorbar;colormap('jet');
subplot(2,2,4);contour(ww,ww,abs(cas2)./abs(dir2),50);colorbar;colormap('jet');