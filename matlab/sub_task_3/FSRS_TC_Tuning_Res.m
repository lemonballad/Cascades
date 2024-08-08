function [ratio,cas,dir] = FSRS_TC_Tuning_Res(E_vib,ovlp,nquanta,...
    parameters_material,parameters_laser)
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

%%
iq=nquanta;
% Speed of light in cm/fs
c=2.998E-5;%
% 2pi*c
r2w=2*pi*c;%0.0001885;%
% Vibrational energy gaps
[col,row]=meshgrid(E_vib,E_vib);
w=row-col;

%% Material parameters
% Electronic dephasing rate
gamma_eg=parameters_material.gamma_eg;
% Vibrational dephasing rate
gamma_vib=parameters_material.gamma_vib;
% Electronic energy gap
weg=parameters_material.weg;

%% Laser parameters
% Frequency pulse actinic
w_ap=parameters_laser.w_ap;
% Frequency pulse raman
w_rp=parameters_laser.w_rp;

% Spectral width of actinic pulse
LAMBDA_ap=parameters_laser.LAMBDA_ap;
% Spectral width of raman pulse
LAMBDA_rp=parameters_laser.LAMBDA_rp;
%

%%
damp=gamma_vib*c;

% Time delay defined from 0 to 10^4 fs (10 ps)
tau_1=0;
tau2=300;
% Frequency signal
w_t=w_rp-parameters_material.wvib;%w_rp-(kj-1)*32;

% Initialize response functions to 0
r1=zeros(iq,1,'double');
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1A=zeros(iq,1,'double');
f1B=f1A;f2A=f1A;f2B=f1A;f3A=f1A;f3B=f1A;f4A=f1A;f4B=f1A;

L_ap_p=1./(w_ap-weg-w+1i*gamma_eg);
L_ap_m=1./(-w_ap+weg-w+1i*gamma_eg);
L_rp_m=1./(-w_rp+weg-w+1i*gamma_eg);
L_t=1./(w_t-weg-w+1i*gamma_eg);
% D=exp((-1i*wD*r2w-damp).*tau_1).*2*LAMBDA_ap./(wD.^2+LAMBDA_ap^2);
D=exp((-1i*w*r2w-damp).*tau_1).*2*LAMBDA_ap./(w.^2+LAMBDA_ap^2);
J=1./(w_t-w_rp-w+1i*(gamma_vib-LAMBDA_rp)).*exp(-LAMBDA_rp*tau2*c);

%% Loop over tau
    % Loop over initial state m => assume m is in the ground state m=1, 0
    % quanta.
    for m=1:1
        % Loop over vibrational states n
        parfor n=1:iq
            % Loop over vibrational states k
            for k=1:iq
                % Loop over vibrational states l
                for l=1:iq
                    % Compute third order cascade auxillary functions
                    % Sequential cascades
                    if n~=k
                        % Auxillary function #4 for third order sequential
                        % cascade weg+dw-.wvib
                        f1A(n)=f1A(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                            *L_t(n,m)*J(n,k)*L_t(n,l);
                        % Auxillary function #2 for third order sequential
                        % cascade
                        f2A(n)=f2A(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                            *L_rp_m(m,n)*J(k,n)*L_t(k,l);
                    end
                    if m~=k
                        % Auxillary function #3 for third order sequential
                        % cascade
                        f3A(n)=f3A(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                            *L_t(n,m)*J(k,m)*L_t(l,m);
                    end
                    if m~=k
                        % Auxillary function #4 for third order sequential
                        % cascade
                        f4A(n)=f4A(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                            *L_rp_m(m,n)*J(m,k)*L_t(l,k);
                    end
                                                
                    % Sequential cascades
                    if k~=m
                        % Auxillary function #1 for third order 
                        % cascade
                        f1B(n)=f1B(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                            *L_ap_p(n,m)*L_t(n,l)*D(n,k);

                        % Auxillary function #2 for third order 
                        % cascade
                        f2B(n)=f2B(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                            *L_ap_m(m,n)*L_t(k,l)*D(k,n);
                    end
                    if k~=m
                        % Auxillary function #3 for third order 
                        % cascade
                        f3B(n)=f3B(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                            *L_ap_p(n,m)*L_t(l,m)*D(k,m);
                        % Auxillary function #4 for third order 
                        % cascade
                        f4B(n)=f4B(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                            *L_ap_m(m,n)*L_t(l,k)*D(m,k);
                    end
                    %
                    
                    % Compute fifth order auxillary functions
                    % Loop over vibrational states u
                    for u=1:iq
                        % Loop over vibrational states v
                        for v=1:iq
                            %
                            if m~=u
                                % Auxillary function #1 for fifth order signal
                               r1(n)=r1(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                    *L_ap_p(n,m)*D(k,m)*L_t(l,m)*J(u,m)*L_t(v,m);
                            end
                            %
                            if k~=u
                                % Auxillary function #2 for fifth order signal
                                r2(n)=r2(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                    *L_ap_p(n,m)*D(k,m)*L_rp_m(k,l)*J(k,u)*L_t(v,u);
                            end
                            %                                                          
                            if m~=u
                                % Auxillary function #3 for fifth order signal
                                r3(n)=r3(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                    *L_ap_m(m,n)*D(m,k)*L_rp_m(m,l)*J(m,u)*L_t(v,u);
                            end
                            %
                            if k~=u
                                % Auxillary function #4 for fifth order signal
                                r4(n)=r4(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                    *L_ap_m(m,n)*D(m,k)*L_t(l,k)*J(u,k)*L_t(v,k);
                            end
                            %
                            if k~=u
                                % Auxillary function #5 for fifth order signal
                                r5(n)=r5(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *L_ap_p(n,m)*D(n,k)*L_rp_m(l,k)*J(u,k)*L_t(u,v);
                            end
                            %
                            if n~=u
                                % Auxillary function #6 for fifth order signal
                                r6(n)=r6(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *L_ap_p(n,m)*D(n,k)*L_t(n,l)*J(n,u)*L_t(n,v);
                            end
                            %
                            if n~=u
                                % Auxillary function #7 for fifth order signal
                                r7(n)=r7(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *L_ap_m(m,n)*D(k,n)*L_rp_m(l,n)*J(u,n)*L_t(u,v);
                            end
                            %
                            if k~=u
                                % Auxillary function #8 for fifth order signal
                                r8(n)=r8(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *L_ap_m(m,n)*D(k,n)*L_t(k,l)*J(k,u)*L_t(k,v);
                            end
                            %
                            if l~=u
                                % Auxillary function #9 for fifth order signal
                                r9(n)=r9(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *L_ap_p(n,m)*D(k,m)*L_t(l,m)*J(l,u)*L_t(l,v);
                                %
                                % Auxillary function #10 for fifth order signal
                                r10(n)=r10(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *L_ap_p(n,m)*D(k,m)*L_rp_m(k,l)*J(u,l)*L_t(u,v);
                                %
                                % Auxillary function #11 for fifth order signal
                                r11(n)=r11(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *L_ap_m(m,n)*D(m,k)*L_t(l,k)*J(l,u)*L_t(l,v);
                                %
                                % Auxillary function #12 for fifth order signal
                                r12(n)=r12(n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *L_ap_m(m,n)*D(m,k)*L_rp_m(m,l)*J(u,l)*L_t(u,v);
                                %
                                % Auxillary function #13 for fifth order signal
                                r13(n)=r13(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *L_ap_p(n,m)*D(n,k)*L_rp_m(l,k)*J(l,u)*L_t(v,u);
                                %
                                % Auxillary function #14 for fifth order signal
                                r14(n)=r14(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *L_ap_p(n,m)*D(n,k)*L_t(n,l)*J(u,l)*L_t(v,l);
                                %
                                % Auxillary function #15 for fifth order signal
                                r15(n)=r15(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *L_ap_m(m,n)*D(k,n)*L_rp_m(l,n)*J(l,u)*L_t(v,u);
                                %
                                % Auxillary function #16 for fifth order signal
                                r16(n)=r16(n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
                                    *L_ap_m(m,n)*D(k,n)*L_t(k,l)*J(u,l)*L_t(v,l);
                            end
                            %
                        end % End loop over vibrational states v
                    end % End loop over vibrational states u
                    
                end % End loop over vibrational state l
            end % End loop over vibrational states k
        end % End loop over vibrational states n
        r1=sum(r1);r2=sum(r2);r3=sum(r3);r4=sum(r4);
        r5=sum(r5);r6=sum(r6);r7=sum(r7);r8=sum(r8);
        r9=sum(r9);r10=sum(r10);r11=sum(r11);r12=sum(r12);
        r13=sum(r13);r14=sum(r14);r15=sum(r15);r16=sum(r16);
        f1A=sum(f1A);f2A=sum(f2A);f3A=sum(f3A);f4A=sum(f4A);
        f1B=sum(f1B);f2B=sum(f2B);f3B=sum(f3B);f4B=sum(f4B);
    end % End loop over initial state

% Compute sequential third order cascade response.
cas=-1i^3*(f1A+f2A+f3A+f4A)...
    .*(f1B+f2B+f3B+f4B);

% COMPUTE Direct fifth order auxillary
direct=(r1+r2+r3+r4)+(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)+(r13+r14+r15+r16);
dir=-direct;

ratio=abs(cas)./abs(dir);

end % End function response2_TC