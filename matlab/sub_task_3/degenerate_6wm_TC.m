function [direct,seq2,e3,tau,wrs] = degenerate_6wm_TC(E_vib,fcall,parameters)
% Compute third order cascade response functions
%   wviball             : vibrational energies for each basis state
%   fcall                 : overlap integrals between basis states.
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
iq=4;
% Speed of light in cm/fs
c=2.998E-5;
% 2pi*c
r2w=2*pi*c;
% Vibrational energy gaps
[col,row]=meshgrid(E_vib,E_vib);
w=row-col;

% Electronic dephasing rate
gamma_eg=parameters.gamma_eg;
% Vibrational dephasing rate
gamma_vib=parameters.gamma_vib;
% Electronic energy gap
weg=parameters.weg;

% Frequency pulse actinic
w_ap=weg+0;
% Frequency pulse raman
w_rp=w_ap;

% Spectral width of pulse a
LAMBDA_ap=210;
% Spectral width of pulse b
LAMBDA_rp=40;

tau2=300*c;
damp=0;

%%
kj=1:100;
% Time delay defined from 0 to 10^4 fs (10 ps)
tau=(kj-1)*100;
% Number of time points in time delay
nt=length(tau);
% Frequency signal
w_t=w_rp-(kj-1)*32;
wrs=w_rp-w_t;

% Initialize response functions to 0
r1=zeros(nt,nt,'double');f1A=zeros(nt,nt,'double');
r2=r1;r3=r1;r4=r1;r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;r13=r1;r14=r1;r15=r1;r16=r1;
f1B=f1A;f2A=f1A;f2B=f1A;f3A=f1A;f3B=f1A;f4A=f1A;f4B=f1A;

% Loop over tau
for ii=1:nt
    % Loop over initial state m => assume m is in the ground state m=1, 0
    % quanta.
    for m=1:1
        % Loop over vibrational states n
        for n=1:iq
            % Loop over vibrational states k
            for k=1:iq
                % Loop over vibrational states l
                for l=1:iq
                    %
                    % DO FWM PARTS OF CASCADES
                    % k3-k4+k5
                    %
                    % damp=(gamma_vib)*c;
                    % if k==m
                    % damp=0;
                    % end
                    if n~=k
                        f1A(ii,:)=f1A(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(k,l)*fcall(n,l)...
                            *1./(w_t-weg-w(n,m)+1i*gamma_eg)...
                            *1./(w_t-w_rp-w(n,k)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                            *1./(w_t-weg-w(n,l)+1i*gamma_eg);
                        %
                        f2A(ii,:)=f2A(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(n,l)*fcall(k,l)...
                            *1./(-w_rp+weg-w(m,n)+1i*gamma_eg)...
                            *1./(w_t-w_rp-w(k,n)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                            *1./(w_t-weg-w(k,l)+1i*gamma_eg);
                    end
                    %
                    if m~=k
                        f3A(ii,:)=f3A(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,k)*fcall(l,m)...
                            *1./(w_t-weg-w(n,m)+1i*gamma_eg)...
                            *1./(w_t-w_rp-w(k,m)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                            *1./(w_t-weg-w(l,m)+1i*gamma_eg);
                    end
                    %
                    if l~=k
                        f4A(ii,:)=f4A(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,m)*fcall(l,k)...
                            *1./(-w_rp+weg-w(m,n)+1i*gamma_eg)...
                            *1./(w_t-w_rp-w(l,k)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                            *1./(w_t-weg-w(m,k)+1i*gamma_eg);
                    end
                    %
                    if k~=m
                        f1B(ii,:)=f1B(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(k,l)*fcall(n,l)...
                            *1./(w_ap-weg-w(n,m)+1i*gamma_eg)...
                            *1./(w_t-weg-w(n,l)+1i*gamma_eg)...
                            .*exp(-1i*w(n,k)*r2w*tau(ii)-damp*tau(ii))*2*LAMBDA_ap/(w(n,k)^2+LAMBDA_ap^2);
                        %
                        f2B(ii,:)=f2B(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(n,l)*fcall(k,l)...
                            *1./(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                            *1./(w_t-weg-w(k,l)+1i*gamma_eg)...
                            .*exp(-1i*w(k,n)*r2w*tau(ii)-damp*tau(ii))*2*LAMBDA_ap/(w(k,n)^2+LAMBDA_ap^2);
                        %
                        f3B(ii,:)=f3B(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,k)*fcall(l,m)...
                            *1./(w_ap-weg-w(n,m)+1i*gamma_eg)...
                            *1./(w_t-weg-w(l,m)+1i*gamma_eg)...
                            .*exp(-1i*w(k,m)*r2w*tau(ii)-damp*tau(ii))*2*LAMBDA_ap/(w(k,m)^2+LAMBDA_ap^2);
                        %
                        f4B(ii,:)=f4B(ii,:)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,m)*fcall(l,k)...
                            *1./(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                            *1./(w_t-weg-w(l,k)+1i*gamma_eg)...
                            .*exp(-1i*w(m,k)*r2w*tau(ii)-damp*tau(ii))*2*LAMBDA_ap/(w(m,k)^2+LAMBDA_ap^2);
                    end
                    %
                    %
                    % Loop over vibrational states u
                    for u=1:iq
                        % Loop over vibrational states v
                        for v=1:iq
                            %
                            if m~=u
                                r1(ii,kj)=r1(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,k)*fcall(l,u)*fcall(v,u)*fcall(v,m)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,m)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,m)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(l,m)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,m)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,m)+1i*gamma_eg);
                            end
                            %
                            if k~=u
                                r2(ii,kj)=r2(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,m)*fcall(l,u)*fcall(v,k)*fcall(v,u)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,m)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,m)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(k,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(k,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,u)+1i*gamma_eg);
                            end
                            %
                            if m~=u
                                r3(ii,kj)=r3(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,k)*fcall(l,u)*fcall(v,m)*fcall(v,u)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(m,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(m,k)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(m,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(m,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,u)+1i*gamma_eg);
                            end
                            %
                            if k~=u
                                r4(ii,kj)=r4(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,m)*fcall(l,u)*fcall(v,u)*fcall(v,k)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(m,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(m,k)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(l,k)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,k)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,k)+1i*gamma_eg);
                            end
                            %
                            if k~=u
                                r5(ii,kj)=r5(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(n,l)*fcall(u,l)*fcall(k,v)*fcall(u,v)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(n,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(n,k)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(l,k)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,k)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(u,v)+1i*gamma_eg);
                            end
                            %
                            if n~=u
                                r6(ii,kj)=r6(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(k,l)*fcall(u,l)*fcall(u,v)*fcall(n,v)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(n,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(n,k)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(n,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(n,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(n,v)+1i*gamma_eg);
                            end
                            %
                            if n~=u
                                r7(ii,kj)=r7(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(k,l)*fcall(u,l)*fcall(n,v)*fcall(u,v)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,n)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,n)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(l,n)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,n)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(u,v)+1i*gamma_eg);
                            end
                            %
                            if k~=u
                                r8(ii,kj)=r8(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(n,l)*fcall(u,l)*fcall(u,v)*fcall(k,v)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,n)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,n)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(k,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(k,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(k,v)+1i*gamma_eg);
                            end
                            %
                            if l~=u
                                r9(ii,kj)=r9(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,k)*fcall(u,m)*fcall(u,v)*fcall(l,v)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,m)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,m)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(l,m)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(l,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(l,v)+1i*gamma_eg);
                                %
                                r10(ii,kj)=r10(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,m)*fcall(u,k)*fcall(l,v)*fcall(u,v)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,m)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,m)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(k,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,l)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(u,v)+1i*gamma_eg);
                                %
                                r11(ii,kj)=r11(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,m)*fcall(u,k)*fcall(u,v)*fcall(l,v)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(m,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(m,k)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(l,k)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(l,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(l,v)+1i*gamma_eg);
                                %
                                r12(ii,kj)=r12(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(n,k)*fcall(l,k)*fcall(u,m)*fcall(l,v)*fcall(u,v)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(m,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(m,k)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(m,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,l)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(u,v)+1i*gamma_eg);
                                %
                                r13(ii,kj)=r13(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(n,l)*fcall(k,u)*fcall(v,l)*fcall(v,u)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(n,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(n,k)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(l,k)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(l,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,u)+1i*gamma_eg);
                                %
                                r14(ii,kj)=r14(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(k,l)*fcall(n,u)*fcall(v,u)*fcall(v,l)...
                                    *1/(w_ap-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp(-1i*w(n,k)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(n,k)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(n,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,l)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,l)+1i*gamma_eg);
                                %
                                r15(ii,kj)=r15(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(k,l)*fcall(n,u)*fcall(v,l)*fcall(v,u)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,n)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,n)^2+LAMBDA_ap^2)...
                                    *1./(-w_rp+weg-w(l,n)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(l,u)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,u)+1i*gamma_eg);
                                %
                                r16(ii,kj)=r16(ii,kj)+boltz_pop(m)*fcall(n,m)*fcall(k,m)*fcall(n,l)*fcall(k,u)*fcall(v,u)*fcall(v,l)...
                                    *1/(-w_ap+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp(-1i*w(k,n)*r2w*tau(ii)-damp*tau(ii) )*2*LAMBDA_ap/(w(k,n)^2+LAMBDA_ap^2)...
                                    *1./(w_t-weg-w(k,l)+1i*gamma_eg)...
                                    *1./(w_t-w_rp-w(u,l)+1i*(gamma_vib-LAMBDA_rp))*exp(-LAMBDA_rp*tau2*c)...
                                    *1./(w_t-weg-w(v,l)+1i*gamma_eg);
                            end
                            %
                            %
                            %
                            %
                        end % End loop over vibrational states v
                    end % End loop over vibrational states u
                end % End loop over vibrational state l
            end % End loop over vibrational states k
        end % End loop over vibrational states n
    end % End loop over initial state

    % Compute sequential third order cascade response.
    seq2(ii,:)=-1i^3*(f1A(ii,:)+f2A(ii,:)+f3A(ii,:)+f4A(ii,:)).*(f1B(ii,:)+f2B(ii,:)+f3B(ii,:)+f4B(ii,:)    );
    % Compute third order response of molecule A
    e3=-1i*(f1A(ii,:)+f2A(ii,:)+f3A(ii,:)+f4A(ii,:));
end % End loop over tau

% COMPUTE SEQUENTIAL CASCADE
direct=(r1+r2+r3+r4)+(r5+r6+r7+r8)+(r9+r10+r11+r12)+(r13+r14+r15+r16);
direct=-direct;

end % End function response2_TC