function [ratio,cascade,direct] = cascade_2dRR_OffRes(E_vib,nquanta,ovlp,...
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
% Vibrational energy gaps Solvent
[wSi,wSf]=meshgrid(wsolv,wsolv);
wS=wSf-wSi;

% Initialize terms of response function to 0
r1=complex(zeros(iq,1,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1_Sol=complex(zeros(iq,1,'double'));
f2_Sol=f1_Sol;f3_Sol=f1_Sol;f4_Sol=f1_Sol;
f1_SolC=complex(zeros(iq,1,'double'));
f2_SolC=f1_SolC;f3_SolC=f1_SolC;f4_SolC=f1_SolC;
f1_Solv=complex(zeros(iq,1,'double'));f2_Solv=f1_Solv;
f1_SolvC=complex(zeros(iq,1,'double'));f2_SolvC=f1_SolvC;

dw=1/nt;
ww=(-1/2:dw:(1/2-dw))/dt/c;
nw=length(ww);
[~,iomega]=min(abs(ww-wvib));
f1PA=complex(zeros(nw,iq,'double'));f2PA=f1PA;f3PA=f1PA;f4PA=f1PA;
f1PB=complex(zeros(nw,iq,'double'));f2PB=f1PB;f3PB=f1PB;f4PB=f1PB;
f1PSA=complex(zeros(nw,iq,'double'));f2PSA=f1PSA;
f1PSB=complex(zeros(nw,iq,'double'));f2PSB=f1PSB;

[W,Wf,Wi]=ndgrid(ww,E_vib,E_vib);W=-W+Wf-Wi;
[WS,WSf,WSi]=ndgrid(ww,wsolv,wsolv);WS=-WS+WSf-WSi;

Lp=1./(w_L-weg-w+1i*gamma_eg);
Lm=1./(-w_L+weg-w+1i*gamma_eg);
Lc1=1./(-wvib-w+1i*gamma_vib);%exp((-1i*wD*r2w-damp).*tau_1D);
Lc2=1./(-wvib-w+1i*gamma_vib);%exp((-1i*wD*r2w-damp).*tau_2D);
Lc0=1./(-wvib-w+1i*gamma_vib);%exp((-1i*wC*r2w-damp).*tau0);
Lc0C=1./(-wvib+w+1i*gamma_vib);%exp((-1i*wC*r2w-damp).*tau0);
LcS=1./(-wvib-wS+1i*gamma_vib);%exp((-1i*wS*r2w-damp).*tau0);
LcSC=1./(-wvib+wS+1i*gamma_vib);%exp((-1i*wS*r2w-damp).*tau0);
LcP=1./(W+1i*gamma_vib);%exp((-1i*wS*r2w-damp).*tau0);
LcPS=1./(WS+1i*gamma_vib);%exp((-1i*wS*r2w-damp).*tau0);
%% Loop over initial state m => assume m is in the ground state m=1, 0
% quanta.
parfor m=1:iq
    % Loop over vibrational states n
    for n=1:iq
        if abs(n-m)==1
            %
            f1_Solv(m)=f1_Solv(m)+boltz_pop_Solv(m)*LcS(n,m);
            f2_Solv(m)=f2_Solv(m)+boltz_pop_Solv(m)*LcS(m,n);
            f1_SolvC(m)=f1_SolvC(m)+boltz_pop_Solv(m)*LcSC(n,m);
            f2_SolvC(m)=f2_SolvC(m)+boltz_pop_Solv(m)*LcSC(m,n);
            f1PSA(:,m)=f1PSA(:,m)+boltz_pop_Solv(m)*LcPS(:,n,m);
            f2PSA(:,m)=f2PSA(:,m)+boltz_pop_Solv(m)*LcPS(:,m,n);
            f1PSB(:,m)=f1PSB(:,m)+boltz_pop_Solv(m)*LcPS(:,n,m)*LcS(n,m);
            f2PSB(:,m)=f2PSB(:,m)+boltz_pop_Solv(m)*LcPS(:,m,n)*LcS(m,n);
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
                    f1_Sol(m)=f1_Sol(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*Lc0(k,m)*Lp(l,m);
                    f2_Sol(m)=f2_Sol(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*Lc0(m,k)*Lp(l,k);
                    f1_SolC(m)=f1_SolC(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*Lc0C(k,m)*Lp(l,m);
                    f2_SolC(m)=f2_SolC(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*Lc0C(m,k)*Lp(l,k);
                    f1PA(:,m)=f1PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP(:,k,m)*Lp(l,m);
                    f2PA(:,m)=f2PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP(:,m,k)*Lp(l,k);
                    f1PB(:,m)=f1PB(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP(:,k,m)*Lc0(m,k)*Lp(l,m);
                    f2PB(:,m)=f2PB(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP(:,m,k)*Lc0C(m,k)*Lp(l,k);
                    %
                end
                %
                if k~=n
                    f3_Sol(m)=f3_Sol(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*Lc0(n,k)*Lp(n,l);
                    f4_Sol(m)=f4_Sol(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*Lc0(k,n)*Lp(k,l);
                    f3_SolC(m)=f3_SolC(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*Lc0C(n,k)*Lp(n,l);
                    f4_SolC(m)=f4_SolC(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*Lc0C(k,n)*Lp(k,l);
                    f3PA(:,m)=f3PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lp(n,l);
                    f4PA(:,m)=f4PA(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lp(k,l);
                    f3PB(:,m)=f3PB(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lc0C(n,k)*Lp(n,l);
                    f4PB(:,m)=f4PB(:,m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lc0C(k,n)*Lp(k,l);
                end
                %
                % Compute fifth order auxillary functions
                % Loop over vibrational states u
                for u=1:iq
                    % Loop over vibrational states v
                    for v=1:iq
                        %
                        if m~=u && m~=k
                            % Auxillary function #1 for fifth order signal
%                             r1(m)=r1(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
%                                 *Lp(n,m)...
%                                 .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                 *Lp(l,m)...
%                                 .*exp((-1i*w(u,m)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,m);
                            r1(m)=r1(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                *Lp(n,m)*Lc1(k,m)*Lp(l,m)*Lc2(u,m)*Lp(v,m);
                            
                        end
                        %
                        if k~=u && m~=k
                            % Auxillary function #2 for fifth order signal
%                             r2(m)=r2(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
%                                 *Lp(n,m)...
%                                 .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                 *Lm(k,l)...
%                                 .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,u);
                            r2(m)=r2(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                *Lp(n,m)*Lc1(k,m)*Lm(k,l)*Lc2(k,u)*Lp(v,u);
                        end
                        %
                        if m~=u && m~=k
                            % Auxillary function #3 for fifth order signal
%                             r3(m)=r3(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
%                                 *Lm(m,n)...
%                                 .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                 *Lm(m,l)...
%                                 .*exp((-1i*w(m,u)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,u);
                            r3(m)=r3(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                *Lm(m,n)...
                                .*Lc1(m,k)...
                                *Lm(m,l)...
                                .*Lc2(m,u)...
                                *Lp(v,u);
                        end
                        %
                        if k~=u && m~=k
                            % Auxillary function #4 for fifth order signal
%                             r4(m)=r4(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
%                                 *Lm(m,n)...
%                                 .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                 *Lp(l,k)...
%                                 .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,k);
                            r4(m)=r4(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                *Lm(m,n)...
                                .*Lc1(m,k)...
                                *Lp(l,k)...
                                .*Lc2(u,k)...
                                *Lp(v,k);
                        end
                        %
                        %
                        if true
                            if k~=u && n~=k
                                % Auxillary function #5 for fifth order signal
%                                 r5(m)=r5(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lm(l,k)...
%                                     .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r5(m)=r5(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(u,k)...
                                    *Lp(u,v);
                            end
                            %
                            if n~=u && n~=k
                                % Auxillary function #6 for fifth order signal
%                                 r6(m)=r6(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lp(n,l)...
%                                     .*exp((-1i*w(n,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(n,v);
                                r6(m)=r6(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(n,u)...
                                    *Lp(n,v);
                                %
                                % Auxillary function #7 for fifth order signal
%                                 r7(m)=r7(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lm(l,n)...
%                                     .*exp((-1i*w(u,n)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r7(m)=r7(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(u,n)...
                                    *Lp(u,v);
                            end
                            %
                            if k~=u && n~=k
                                % Auxillary function #8 for fifth order signal
%                                 r8(m)=r8(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lp(k,l)...
%                                     .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(k,v);
                                r8(m)=r8(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lp(k,l)...
                                    .*Lc2(k,u)...
                                    *Lp(k,v);
                            end
                            %
                            if l~=u && m~=k
                                % Auxillary function #9 for fifth order signal
%                                 r9(m)=r9(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                     *Lp(l,m)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(l,v);
                                r9(m)=r9(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *Lp(n,m)...
                                    .*Lc1(k,m)...
                                    *Lp(l,m)...
                                    .*Lc2(l,u)...
                                    *Lp(l,v);
                                %
                                % Auxillary function #10 for fifth order signal
%                                 r10(m)=r10(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                     *Lm(k,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r10(m)=r10(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    .*Lc1(k,m)...
                                    *Lm(k,l)...
                                    .*Lc2(u,l)...
                                    *Lp(u,v);
                                %
                                % Auxillary function #11 for fifth order signal
%                                 r11(m)=r11(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                     *Lp(l,k)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(l,v);
                                r11(m)=r11(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *Lm(m,n)...
                                    .*Lc1(m,k)...
                                    *Lp(l,k)...
                                    .*Lc2(l,u)...
                                    *Lp(l,v);
                                %
                                % Auxillary function #12 for fifth order signal
%                                 r12(m)=r12(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                     *Lm(m,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r12(m)=r12(m)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    .*Lc1(m,k)...
                                    *Lm(m,l)...
                                    .*Lc2(u,l)...
                                    *Lp(u,v);
                            end
                            if l~=u && n~=k
                                %
                                % Auxillary function #13 for fifth order signal
%                                 r13(m)=r13(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lm(l,k)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,u);
                                r13(m)=r13(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(l,u)...
                                    *Lp(v,u);
                                %
                                % Auxillary function #14 for fifth order signal
%                                 r14(m)=r14(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lp(n,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,l);
                                r14(m)=r14(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(u,l)...
                                    *Lp(v,l);
                                %
                                % Auxillary function #15 for fifth order signal
%                                 r15(m)=r15(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lm(l,n)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,u);
                                r15(m)=r15(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(l,u)...
                                    *Lp(v,u);
                                %
                                % Auxillary function #16 for fifth order signal
%                                 r16(m)=r16(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lp(k,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,l);
                                r16(m)=r16(m)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lp(k,l)...
                                    .*Lc2(u,l)...
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
f1_Sol=sum(f1_Sol);f2_Sol=sum(f2_Sol);f3_Sol=sum(f3_Sol);f4_Sol=sum(f4_Sol);
f1_SolC=sum(f1_SolC);f2_SolC=sum(f2_SolC);f3_SolC=sum(f3_SolC);f4_SolC=sum(f4_SolC);
f1PA=sum(f1PA,2);f2PA=sum(f2PA,2);f3PA=sum(f3PA,2);f4PA=sum(f4PA,2);
f1PB=sum(f1PB,2);f2PB=sum(f2PB,2);f3PB=sum(f3PB,2);f4PB=sum(f4PB,2);
f1_Solv=sum(f1_Solv);f2_Solv=sum(f2_Solv);
f1_SolvC=sum(f1_SolvC);f2_SolvC=sum(f2_SolvC);
f1PSA=sum(f1PSA,2);f2PSA=sum(f2PSA,2);
f1PSB=sum(f1PSB,2);f2PSB=sum(f2PSB,2);
r1=sum(r1);r2=sum(r2);r3=sum(r3);r4=sum(r4);
r5=sum(r5);r6=sum(r6);r7=sum(r7);r8=sum(r8);
r9=sum(r9);r10=sum(r10);r11=sum(r11);r12=sum(r12);
r13=sum(r13);r14=sum(r14);r15=sum(r15);r16=sum(r16);

%
direct=((r1+r2+r3+r4)+(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)+(r13+r14+r15+r16));

% COMPUTE 2D CASCADES
seq1=-(f1_Solv+f2_Solv+f1_SolvC+f2_SolvC)...
    *(f1_Sol+f2_Sol+f3_Sol+f4_Sol);
seq2=-(f1_Sol+f2_Sol+f3_Sol+f4_Sol+f1_SolC+f2_SolC+f3_SolC+f4_SolC)...
    *(f1_Solv+f2_Solv);
wr=1/nw;
par1=-2*(conv((f1PA+f2PA+f3PA+f4PA),-1i*(f1PSB+f2PSB),'same'))*wr;
par2=-2*(conv((f1PSA+f2PSA),-1i*(f1PB+f2PB+f3PB+f4PB),'same'))*wr;

seq=seq1+seq2;
par=par1+par2;
cascade=seq+par(iomega);
ratio=abs(cascade)/abs(direct);



