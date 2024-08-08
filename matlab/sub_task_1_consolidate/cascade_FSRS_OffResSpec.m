function [cas2,dir2] = cascade_FSRS_OffResSpec(wviball,nquanta,ovlp,...
    parameters_laser,parameters_material)
% Compute third order cascade response functions
%   wviball             : vibrational energies for each basis state
%   ovlp                  : overlap integrals between basis states.
%   nquanta           : number of quanta
%   parameters      : Material parameters

%%
iq=nquanta;
% Speed of light in cm/fs
c=2.998E-5;%
r2w=2*pi*c;
% Vibrational energy gaps
[wi,wf]=meshgrid(wviball(1:nquanta+1),wviball(1:nquanta+1));
w=wf-wi;

%% Material parameters
% Electronic dephasing rate
gamma_eg=parameters_material.gamma_eg;
% Vibrational dephasing rate
gamma_vib=parameters_material.gamma_vib;
% Electronic energy gap
weg=parameters_material.weg;
% Vibrational frequencies Solvent
wSolv=parameters_material.wsolv;
% Vibrational energy levels Solvent
wsolv=(0:nquanta+1/2)*wSolv;

%% COMPUTE BOLTZMANN POPULATIONS
kT=200; % kT in cm^-1
% Compute array of
boltz_factor=exp(-wviball/kT);
boltz_factor_Solv=exp(-wsolv/kT);
% Compute partition function
partition_func=sum(boltz_factor);
partition_func_Solv=sum(boltz_factor_Solv);
% Compute boltzman populations
boltz_pop=boltz_factor/partition_func;
boltz_pop_Solv=boltz_factor_Solv/partition_func_Solv;

%% Laser parameters
% Frequency pulse actinic
w_ap=parameters_laser.w_ap;
% Frequency pulse raman
w_rp=parameters_laser.w_rp;

% Spectral width of actinic pulse
LAMBDA_ap=parameters_laser.LAMBDA_ap;
% Spectral width of raman pulse
LAMBDA_rp=parameters_laser.LAMBDA_rp;
% Time step
dt=parameters_laser.dt;
% Number of time steps
nt=parameters_laser.nt;

%%
damp=gamma_vib*r2w;

% Time delay
tau_1=0:dt:(nt-1)*dt;
tau_2=0;
% Number of frequency points in spectrum
nw=nt;
% Frequency step
dw=1/nw;
% Frequency axis
ww=(-1/2:dw:(1/2-dw))/dt/c;ww=ww-ww(1);

% Initialize response functions to 0
r1=complex(zeros(nw,nt,iq,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1A=complex(zeros(nw,nt,iq,'double'));
f1B=f1A;f2A=f1A;f2B=f1A;f3A=f1A;f3B=f1A;f4A=f1A;f4B=f1A;
f1SA=complex(zeros(nw,nt,iq,'double'));
f1SB=f1SA;f2SA=f1SA;f2SB=f1SA;

[W,Tau_1,wDf,wDi]=ndgrid(ww,tau_1,wviball(1:nquanta+1),wviball(1:nquanta+1));
wD=wDf-wDi;
DAMP=damp*ones(size(wD));DAMP(wD==0)=0;
W_t=w_rp-W;
[~,~,wSf,wSi]=ndgrid(ww,tau_1,wsolv,wsolv);
wS=wSf-wSi;
L_ap_p=1./(w_ap-weg-w+1i*gamma_eg);
L_ap_m=1./(-w_ap+weg-w+1i*gamma_eg);
L_rp_m=1./(-w_rp+weg-w+1i*gamma_eg);
L_t=1./(W_t-weg-wD+1i*gamma_eg);
D=exp((-1i*wD*r2w-damp).*Tau_1).*(2*LAMBDA_ap)./(wD.^2+LAMBDA_ap^2);
D0=exp((-1i*wD*r2w-DAMP).*Tau_1).*(2*LAMBDA_ap)./(wD.^2+LAMBDA_ap^2);
J=1./(W_t-w_rp-wD+1i*(gamma_vib-LAMBDA_rp)).*exp((1i*(w_rp-W_t)-LAMBDA_rp)*tau_2*r2w);
DS=exp((-1i*wS*r2w-damp).*Tau_1).*(2*LAMBDA_ap)./(wS.^2+LAMBDA_ap^2);
JS=1./(W_t-w_rp-wS+1i*(gamma_vib-LAMBDA_rp)).*exp((1i*(w_rp-W_t)-LAMBDA_rp)*tau_2*r2w);
D(wD==0)=D0(wD==0);
J(wD==0)=0;JS(wS==0)=0;

%% Loop over tau
% Loop over initial state m => assume m is in the ground state m=1, 0
% quanta.
for m=1:1
    % Loop over vibrational states n
    parfor n=1:iq
        % Loop over vibrational states k
        if abs(m-n)==1
            % Auxillary function #1 for third order cascade
            f1SA(:,:,n)=f1SA(:,:,n)+boltz_pop_Solv(m)*JS(:,:,n,m);
            % Auxillary function #2 for third order cascade
            f2SA(:,:,n)=f2SA(:,:,n)+boltz_pop_Solv(m)*JS(:,:,m,n);
            % Auxillary function #1 for third order cascade
            f1SB(:,:,n)=f1SB(:,:,n)+boltz_pop_Solv(m)*DS(:,:,n,m);
            % Auxillary function #2 for third order cascade
            f2SB(:,:,n)=f2SB(:,:,n)+boltz_pop_Solv(m)*DS(:,:,m,n);
        end
        
        for k=1:iq
            % Loop over vibrational states l
            for l=1:iq
                % Compute third order cascade auxillary functions
                % Auxillary function #4 for third order cascade
                f1A(:,:,n)=f1A(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                    .*L_t(:,:,n,m).*J(:,:,n,k).*L_t(:,:,n,l);
                % Auxillary function #2 for third order cascade
                f2A(:,:,n)=f2A(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                    *L_rp_m(m,n).*J(:,:,k,n).*L_t(:,:,k,l);
                % Auxillary function #3 for third order cascade
                f3A(:,:,n)=f3A(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                    .*L_t(:,:,n,m).*J(:,:,k,m).*L_t(:,:,l,m);
                % Auxillary function #4 for third order cascade
                f4A(:,:,n)=f4A(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                    *L_rp_m(m,n).*J(:,:,m,k).*L_t(:,:,l,k);
                % Auxillary function #1 for third order cascade
                f1B(:,:,n)=f1B(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                    *L_ap_p(n,m)*L_t(:,:,n,l).*D(:,:,n,k);
                % Auxillary function #2 for third order cascade
                f2B(:,:,n)=f2B(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                    *L_ap_m(m,n)*L_t(:,:,k,l).*D(:,:,k,n);
                % Auxillary function #3 for third order cascade
                f3B(:,:,n)=f3B(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                    *L_ap_p(n,m).*L_t(:,:,l,m).*D(:,:,k,m);
                % Auxillary function #4 for third order cascade
                f4B(:,:,n)=f4B(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                    *L_ap_m(m,n).*L_t(:,:,l,k).*D(:,:,m,k);
                % Compute fifth order auxillary functions
                % Loop over vibrational states u
                for u=1:iq
                    % Loop over vibrational states v
                    for v=1:iq
                        % Auxillary function #1 for fifth order signal
                        r1(:,:,n)=r1(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                            *L_ap_p(n,m).*D(:,:,k,m).*L_t(:,:,l,m).*J(:,:,u,m).*L_t(:,:,v,m);
                        % Auxillary function #2 for fifth order signal
                        r2(:,:,n)=r2(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                            *L_ap_p(n,m).*D(:,:,k,m)*L_rp_m(k,l).*J(:,:,k,u).*L_t(:,:,v,u);
                        % Auxillary function #3 for fifth order signal
                        r3(:,:,n)=r3(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                            *L_ap_m(m,n).*D(:,:,m,k)*L_rp_m(m,l).*J(:,:,m,u).*L_t(:,:,v,u);
                        % Auxillary function #4 for fifth order signal
                        r4(:,:,n)=r4(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                            *L_ap_m(m,n).*D(:,:,m,k).*L_t(:,:,l,k).*J(:,:,u,k).*L_t(:,:,v,k);
                        % Auxillary function #5 for fifth order signal
                        r5(:,:,n)=r5(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                            *L_ap_p(n,m).*D(:,:,n,k)*L_rp_m(l,k).*J(:,:,u,k).*L_t(:,:,u,v);
                        % Auxillary function #6 for fifth order signal
                        r6(:,:,n)=r6(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                            *L_ap_p(n,m).*D(:,:,n,k).*L_t(:,:,n,l).*J(:,:,n,u).*L_t(:,:,n,v);
                        % Auxillary function #7 for fifth order signal
                        r7(:,:,n)=r7(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                            *L_ap_m(m,n).*D(:,:,k,n)*L_rp_m(l,n).*J(:,:,u,n) .*L_t(:,:,u,v);
                        % Auxillary function #8 for fifth order signal
                        r8(:,:,n)=r8(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                            *L_ap_m(m,n).*D(:,:,k,n).*L_t(:,:,k,l).*J(:,:,k,u).*L_t(:,:,k,v);
                        % Auxillary function #9 for fifth order signal
                        r9(:,:,n)=r9(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                            *L_ap_p(n,m).*D(:,:,k,m).*L_t(:,:,l,m).*J(:,:,l,u).*L_t(:,:,l,v);
                        % Auxillary function #10 for fifth order signal
                        r10(:,:,n)=r10(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                            *L_ap_p(n,m).*D(:,:,k,m)*L_rp_m(k,l).*J(:,:,u,l).*L_t(:,:,u,v);
                        % Auxillary function #11 for fifth order signal
                        r11(:,:,n)=r11(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                            *L_ap_m(m,n).*D(:,:,m,k).*L_t(:,:,l,k).*J(:,:,l,u).*L_t(:,:,l,v);
                        % Auxillary function #12 for fifth order signal
                        r12(:,:,n)=r12(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                            *L_ap_m(m,n).*D(:,:,m,k)*L_rp_m(m,l).*J(:,:,u,l).*L_t(:,:,u,v);
                        % Auxillary function #13 for fifth order signal
                        r13(:,:,n)=r13(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                            *L_ap_p(n,m).*D(:,:,n,k)*L_rp_m(l,k).*J(:,:,l,u).*L_t(:,:,v,u);
                        % Auxillary function #14 for fifth order signal
                        r14(:,:,n)=r14(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                            *L_ap_p(n,m).*D(:,:,n,k).*L_t(:,:,n,l).*J(:,:,u,l).*L_t(:,:,v,l);
                        % Auxillary function #15 for fifth order signal
                        r15(:,:,n)=r15(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                            *L_ap_m(m,n).*D(:,:,k,n)*L_rp_m(l,n).*J(:,:,l,u).*L_t(:,:,v,u);
                        % Auxillary function #16 for fifth order signal
                        r16(:,:,n)=r16(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
                            *L_ap_m(m,n).*D(:,:,k,n).*L_t(:,:,k,l).*J(:,:,u,l).*L_t(:,:,v,l);
                    end % End loop over vibrational states v
                end % End loop over vibrational states u
            end % End loop over vibrational state l
        end % End loop over vibrational states k
    end % End loop over vibrational states n
end % End loop over initial state

r1=sum(r1,3);r2=sum(r2,3);r3=sum(r3,3);r4=sum(r4,3);
r5=sum(r5,3);r6=sum(r6,3);r7=sum(r7,3);r8=sum(r8,3);
r9=sum(r9,3);r10=sum(r10,3);r11=sum(r11,3);r12=sum(r12,3);
r13=sum(r13,3);r14=sum(r14,3);r15=sum(r15,3);r16=sum(r16,3);
f1A=sum(f1A,3);f2A=sum(f2A,3);f3A=sum(f3A,3);f4A=sum(f4A,3);
f1B=sum(f1B,3);f2B=sum(f2B,3);f3B=sum(f3B,3);f4B=sum(f4B,3);
f1SA=sum(f1SA,3);f2SA=sum(f2SA,3);
f1SB=sum(f1SB,3);f2SB=sum(f2SB,3);

% Compute sequential third order cascade response.
seq1=(1i)^4*(f1A+f2A+f3A+f4A)...
    .*(f1SB+f2SB);
seq2=(1i)^4*(f1SA+f2SA)...
    .*(f1B+f2B+f3B+f4B);
cas2=seq1+seq2;

% COMPUTE Direct fifth order auxillary
direct=(1i)^4*((r1+r2+r3+r4)...
    +(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)...
    +(r13+r14+r15+r16));
dir2=direct;

% figure;plot(ww,abs(dir2(:,77)),'k-');
% figure;mesh(tau1,ww,abs(dir2));
% figure;
% subplot(2,2,1);contour(tau_1,ww,8.6298e-16*3e10*abs(cas2),50);colorbar;colormap('jet');title('|E_{Cascade}|');
% subplot(2,2,2);contour(tau_1,ww,abs(dir2),50);colorbar;colormap('jet');title('|E_{Direct}|');
% subplot(2,2,3);contour(tau_1,ww,abs(8.6298e-16*3e10*cas2+dir2),50);colorbar;colormap('jet');title('|E_{Cascade}+E_{Direct}|');
% subplot(2,2,4);contour(tau_1,ww,abs(8.6298e-16*3e10*cas2)./abs(dir2),50);colorbar;colormap('jet');title('|E_{Cascade}|/|E_{Direct}|');
end % End function response2_TC