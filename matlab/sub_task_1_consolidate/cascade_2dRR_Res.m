function [ratio,cas,dir] = cascade_2dRR_Res(E_vib,nquanta,ovlp,...
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
%   direct2d          : Direct fif true%th order signal

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

%
% Initialize terms of response function to 0
r1=complex(zeros(iq,1,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1=complex(zeros(iq,1,'double'));
f2=f1;f3=f1;f4=f1;
f1C=complex(zeros(iq,1,'double'));
f2C=f1C;f3C=f1C;f4C=f1C;

dw=1/nt;
ww=(-1/2:dw:(1/2-dw))/dt/c;
nw=length(ww);
f1PA=complex(zeros(nw,iq,'double'));
f2PA=f1PA;f3PA=f1PA;f4PA=f1PA;
f1P=complex(zeros(nw,iq,'double'));
f2P=f1P;f3P=f1P;f4P=f1P;
[W,Wf,Wi]=ndgrid(ww,E_vib,E_vib);
W=-W+Wf-Wi;
Lp=1./(w_L-weg-w+1i*gamma_eg);
Lm=1./(-w_L+weg-w+1i*gamma_eg);
Lc1=1./(-wvib-w+1i*gamma_vib);%2*gamma_vib./((wvib+w).^2+gamma_vib^2);%1./(wvib+w+1i*gamma_vib);%exp((-1i*wD*r2w-damp).*tau_1D);
Lc2=1./(-wvib-w+1i*gamma_vib);%2*gamma_vib./((wvib+w).^2+gamma_vib^2);%1./(wvib+w+1i*gamma_vib);%exp((-1i*wD*r2w-damp).*tau_2D);
Lc0=1./(-wvib-w+1i*gamma_vib);%2*gamma_vib./((wvib+w).^2+gamma_vib^2);%1./(wvib+w+1i*gamma_vib);%exp((-1i*wC*r2w-damp).*tau0);
Lc0C=1./(wvib+w+1i*gamma_vib);%2*gamma_vib./((wvib-w).^2+gamma_vib^2);%1./(wvib-w+1i*gamma_vib);%exp((-1i*wC*r2w-damp).*tau0);
LcP=1./(-W+1i*gamma_vib);%2*gamma_vib./((wvib-w).^2+gamma_vib^2);%1./(wvib-w+1i*gamma_vib);%exp((-1i*wC*r2w-damp).*tau0);

Lc1(w==0)=0;%1/wvib;
Lc2(w==0)=0;%1/wvib;
Lc0(w==0)=0;%1/wvib;
Lc0C(w==0)=0;%-1/wvib;
LcP(w==0)=0;%1/wvib;
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
                if true% k~=m
                    f1(n)=f1(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*Lc0(k,m)*Lp(l,m);
                    f2(n)=f2(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*Lc0(m,k)*Lp(l,k);
                    f1C(n)=f1C(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*Lc0C(k,m)*Lp(l,m);
                    f2C(n)=f2C(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*Lc0C(m,k)*Lp(l,k);
                    f1PA(:,n)=f1PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP(:,k,m)*Lp(l,m);
                    f2PA(:,n)=f2PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP(:,m,k)*Lp(l,k);
                    f1P(:,n)=f1P(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *Lp(n,m)*LcP(:,k,m)*Lc0(k,m)*Lp(l,m);
                    f2P(:,n)=f2P(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *Lm(m,n)*LcP(:,m,k)*Lc0(m,k)*Lp(l,k);
                end
                %
                if true% k~=n
                    f3(n)=f3(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*Lc0(n,k)*Lp(n,l);
                    f4(n)=f4(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*Lc0(k,n)*Lp(k,l);
                    f3C(n)=f3C(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*Lc0C(n,k)*Lp(n,l);
                    f4C(n)=f4C(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*Lc0C(k,n)*Lp(k,l);
                    f3PA(:,n)=f3PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lp(n,l);
                    f4PA(:,n)=f4PA(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lp(k,l);
                    f3P(:,n)=f3P(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *Lp(n,m)*LcP(:,n,k)*Lc0(n,k)*Lp(n,l);
                    f4P(:,n)=f4P(:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *Lm(m,n)*LcP(:,k,n)*Lc0(k,n)*Lp(k,l);
                end
                %
                % Compute fif true%th order auxillary functions
                % Loop over vibrational states u
                for u=1:iq
                    % Loop over vibrational states v
                    for v=1:iq
                        %
                        if true% m~=u && m~=k
                            % Auxillary function #1 for fif true%th order signal
%                             r1(n)=r1(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
%                                 *Lp(n,m)...
%                                 .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                 *Lp(l,m)...
%                                 .*exp((-1i*w(u,m)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,m);
                            r1(n)=r1(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                *Lp(n,m)...
                                .*Lc1(k,m)...
                                *Lp(l,m)...
                                .*Lc2(u,m)...
                                *Lp(v,m);
                            
                        end
                        %
                        if true% k~=u && m~=k
                            % Auxillary function #2 for fif true%th order signal
%                             r2(n)=r2(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
%                                 *Lp(n,m)...
%                                 .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                 *Lm(k,l)...
%                                 .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,u);
                            r2(n)=r2(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                *Lp(n,m)...
                                .*Lc1(k,m)...
                                *Lm(k,l)...
                                .*Lc2(k,u)...
                                *Lp(v,u);
                        end
                        %
                        if true% m~=u && m~=k
                            % Auxillary function #3 for fif true%th order signal
%                             r3(n)=r3(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
%                                 *Lm(m,n)...
%                                 .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                 *Lm(m,l)...
%                                 .*exp((-1i*w(m,u)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,u);
                            r3(n)=r3(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                *Lm(m,n)...
                                .*Lc1(m,k)...
                                *Lm(m,l)...
                                .*Lc2(m,u)...
                                *Lp(v,u);
                        end
                        %
                        if true% k~=u && m~=k
                            % Auxillary function #4 for fif true%th order signal
%                             r4(n)=r4(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
%                                 *Lm(m,n)...
%                                 .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                 *Lp(l,k)...
%                                 .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
%                                 *Lp(v,k);
                            r4(n)=r4(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                *Lm(m,n)...
                                .*Lc1(m,k)...
                                *Lp(l,k)...
                                .*Lc2(u,k)...
                                *Lp(v,k);
                        end
                        %
                        %
                        if true% true
                            if true% k~=u && n~=k
                                % Auxillary function #5 for fif true%th order signal
%                                 r5(n)=r5(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lm(l,k)...
%                                     .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r5(n)=r5(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(u,k)...
                                    *Lp(u,v);
                            end
                            %
                            if true% n~=u && n~=k
                                % Auxillary function #6 for fif true%th order signal
%                                 r6(n)=r6(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lp(n,l)...
%                                     .*exp((-1i*w(n,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(n,v);
                                r6(n)=r6(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(n,u)...
                                    *Lp(n,v);
                                %
                                % Auxillary function #7 for fif true%th order signal
%                                 r7(n)=r7(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lm(l,n)...
%                                     .*exp((-1i*w(u,n)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r7(n)=r7(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(u,n)...
                                    *Lp(u,v);
                            end
                            %
                            if true% k~=u && n~=k
                                % Auxillary function #8 for fif true%th order signal
%                                 r8(n)=r8(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lp(k,l)...
%                                     .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(k,v);
                                r8(n)=r8(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lp(k,l)...
                                    .*Lc2(k,u)...
                                    *Lp(k,v);
                            end
                            %
                            if true% l~=u && m~=k
                                % Auxillary function #9 for fif true%th order signal
%                                 r9(n)=r9(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                     *Lp(l,m)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(l,v);
                                r9(n)=r9(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *Lp(n,m)...
                                    .*Lc1(k,m)...
                                    *Lp(l,m)...
                                    .*Lc2(l,u)...
                                    *Lp(l,v);
                                %
                                % Auxillary function #10 for fif true%th order signal
%                                 r10(n)=r10(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
%                                     *Lm(k,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r10(n)=r10(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *Lp(n,m)...
                                    .*Lc1(k,m)...
                                    *Lm(k,l)...
                                    .*Lc2(u,l)...
                                    *Lp(u,v);
                                %
                                % Auxillary function #11 for fif true%th order signal
%                                 r11(n)=r11(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                     *Lp(l,k)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(l,v);
                                r11(n)=r11(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *Lm(m,n)...
                                    .*Lc1(m,k)...
                                    *Lp(l,k)...
                                    .*Lc2(l,u)...
                                    *Lp(l,v);
                                %
                                % Auxillary function #12 for fif true%th order signal
%                                 r12(n)=r12(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
%                                     *Lm(m,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(u,v);
                                r12(n)=r12(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *Lm(m,n)...
                                    .*Lc1(m,k)...
                                    *Lm(m,l)...
                                    .*Lc2(u,l)...
                                    *Lp(u,v);
                            end
                            if true% l~=u && n~=k
                                %
                                % Auxillary function #13 for fif true%th order signal
%                                 r13(n)=r13(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lm(l,k)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,u);
                                r13(n)=r13(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lm(l,k)...
                                    .*Lc2(l,u)...
                                    *Lp(v,u);
                                %
                                % Auxillary function #14 for fif true%th order signal
%                                 r14(n)=r14(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
%                                     *Lp(n,m)...
%                                     .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
%                                     *Lp(n,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,l);
                                r14(n)=r14(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *Lp(n,m)...
                                    .*Lc1(n,k)...
                                    *Lp(n,l)...
                                    .*Lc2(u,l)...
                                    *Lp(v,l);
                                %
                                % Auxillary function #15 for fif true%th order signal
%                                 r15(n)=r15(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lm(l,n)...
%                                     .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,u);
                                r15(n)=r15(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *Lm(m,n)...
                                    .*Lc1(k,n)...
                                    *Lm(l,n)...
                                    .*Lc2(l,u)...
                                    *Lp(v,u);
                                %
                                % Auxillary function #16 for fif true%th order signal
%                                 r16(n)=r16(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
%                                     *Lm(m,n)...
%                                     .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
%                                     *Lp(k,l)...
%                                     .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
%                                     *Lp(v,l);
                                r16(n)=r16(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
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
    f1=sum(f1);f2=sum(f2);f3=sum(f3);f4=sum(f4);
    f1C=sum(f1C);f2C=sum(f2C);f3C=sum(f3C);f4C=sum(f4C);
    f1PA=sum(f1PA,2);f2PA=sum(f2PA,2);f3PA=sum(f3PA,2);f4PA=sum(f4PA,2);
    f1P=sum(f1P,2);f2P=sum(f2P,2);f3P=sum(f3P,2);f4P=sum(f4P,2);
    r1=sum(r1);r2=sum(r2);r3=sum(r3);r4=sum(r4);
    r5=sum(r5);r6=sum(r6);r7=sum(r7);r8=sum(r8);
    r9=sum(r9);r10=sum(r10);r11=sum(r11);r12=sum(r12);
    r13=sum(r13);r14=sum(r14);r15=sum(r15);r16=sum(r16);

%
direct=((r1+r2+r3+r4)+(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)+(r13+r14+r15+r16));%0.09*
%
% COMPUTE 2D CASCADES
%
seq1=(f1.^2+f2.^2+f3.^2+f4.^2)...
    +2*(f1.*f2+f1.*f3+f1.*f4+f2.*f3+f2.*f4+f3.*f4);
seq2=(f1C.*f1+f2C.*f2+f3C.*f3+f4C.*f4)...
    +1*(f1C.*f2+f1C.*f3+f1C.*f4+f2C.*f1+f2C.*f3+f2C.*f4...
    +f3C.*f1+f3C.*f2+f3C.*f4+f4C.*f1+f4C.*f2+f4C.*f3);
wr=1/nw;
par1=(conv((f1PA+f2PA+f3PA+f4PA),(f1P+f2P+f3P+f4P),'same'))*wr;
par2=par1;

seq=seq1+seq2;

par=par1+par2;
[~,iomega]=min(abs(ww-wvib));
cascade=seq+par(iomega);

cas=abs(cascade);
dir=abs(direct);
ratio=cas/dir;
