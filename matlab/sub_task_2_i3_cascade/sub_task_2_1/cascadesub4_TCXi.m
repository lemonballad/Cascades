function [ratio,cas,dir] = cascadesub4_TC3(E_vib,nquanta,ovlp,...
    disp,gamma_eg,gamma_vib,weg,dt,nt,w_L,wvib)
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
nt=nt*2;

%% Material parameters
% Electronic dephasing rate
% Vibrational dephasing rate
% Electronic energy gap
% Vibrational frequencies

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
kj=1:nt;
tau=(kj-1)*dt;

% Initialize terms of response function to 0
r1=complex(zeros(nt,nt,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1=complex(zeros(1,nt,'double'));
f2=f1;f3=f1;f4=f1;


[tau_1,tau_2]=meshgrid(tau);
%% Loop over initial state m => assume m is in the ground state m=1, 0
% quanta.
for m=1:1
    % Loop over vibrational states n
    for n=1:iq
        % Loop over vibrational states k
        for k=1:iq
            % Loop over vibrational states l
            for l=1:iq
                %
                % NOTE OVERLAP INEGRALS ARE STORED WITH EXCITED STATE INDEX FIRST
                %
                if true%k~=m
                    f1=f1+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                        .*exp((-1i*w(k,m)*r2w-gamma_vib*c)*tau)...
                        *1/(w_L-weg-w(l,m)+1i*gamma_eg);
                    %
                    f2=f2+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                        .*exp((-1i*w(m,k)*r2w-gamma_vib*c)*tau)...
                        *1/(w_L-weg-w(l,k)+1i*gamma_eg);
                end
                %
                if false%k~=n
                    f3=f3+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *1./(w_L-weg-w(n,m)+1i*gamma_eg)...
                        .*exp((-1i*w(n,k)*r2w-gamma_vib*c)*tau)...
                        *1./(w_L-weg-w(n,l)+1i*gamma_eg);
                    %
                    f4=f4+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
                        *1./(-w_L+weg-w(m,n)+1i*gamma_eg)...
                        .*exp((-1i*w(k,n)*r2w-gamma_vib*c)*tau)...
                        *1./(w_L-weg-w(k,l)+1i*gamma_eg);
                end
                %
                % Compute fifth order auxillary functions
                % Loop over vibrational states u
                for u=1:iq
                    % Loop over vibrational states v
                    for v=1:iq
                        %
                        if true%m~=u && m~=k
                            % Auxillary function #1 for fifth order signal
                            r1=r1+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                *1./(w_L-weg-w(l,m)+1i*gamma_eg)...
                                .*exp((-1i*w(u,m)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,m)+1i*gamma_eg);
                            
                        end
                        %
                        if true%k~=u && m~=k
                            % Auxillary function #2 for fifth order signal
                            r2=r2+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                *1./(-w_L+weg-w(k,l)+1i*gamma_eg)...
                                .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                        end
                        %
                        if true%m~=u && m~=k
                            % Auxillary function #3 for fifth order signal
                            r3=r3+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                *1./(-w_L+weg-w(m,l)+1i*gamma_eg)...
                                .*exp((-1i*w(m,u)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                        end
                        %
                        if true%k~=u && m~=k
                            % Auxillary function #4 for fifth order signal
                            r4=r4+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                *1./(w_L-weg-w(l,k)+1i*gamma_eg)...
                                .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,k)+1i*gamma_eg);
                        end
                        %
                        %
                        if false
                            if k~=u && n~=k
                                % Auxillary function #5 for fifth order signal
                                r5=r5+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,k)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                            end
                            %
                            if n~=u && n~=k
                                % Auxillary function #6 for fifth order signal
                                r6=r6+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(n,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(n,v)+1i*gamma_eg);
                                %
                                % Auxillary function #7 for fifth order signal
                                r7=r7+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,n)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                            end
                            %
                            if k~=u && n~=k
                                % Auxillary function #8 for fifth order signal
                                r8=r8+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(k,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(k,v)+1i*gamma_eg);
                            end
                            %
                            if l~=u && m~=k
                                % Auxillary function #9 for fifth order signal
                                r9=r9+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(l,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(l,v)+1i*gamma_eg);
                                %
                                % Auxillary function #10 for fifth order signal
                                r10=r10+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(k,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                                %
                                % Auxillary function #11 for fifth order signal
                                r11=r11+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(l,k)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(l,v)+1i*gamma_eg);
                                %
                                % Auxillary function #12 for fifth order signal
                                r12=r12+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(m,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                            end
                            if l~=u && n~=k
                                %
                                % Auxillary function #13 for fifth order signal
                                r13=r13+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,k)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                                %
                                % Auxillary function #14 for fifth order signal
                                r14=r14+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(n,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,l)+1i*gamma_eg);
                                %
                                % Auxillary function #15 for fifth order signal
                                r15=r15+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                                %
                                % Auxillary function #16 for fifth order signal
                                r16=r16+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(k,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,l)+1i*gamma_eg);
                            end
                        end
                        %
                        
                        
                        
                    end % End loop over vibrational states v
                end % End loop over vibrational states u
                
            end % End loop over vibrational state l
        end % End loop over vibrational states k
    end % End loop over vibrational states n
end % End loop over initial state

%
% direct=-0.09*(r1+r2+r3+r4);
direct=-((r1+r2+r3+r4));%+(r5+r6+r7+r8)...    +(r9+r10+r11+r12)+(r13+r14+r15+r16));%0.09*
% direct=direct-mean(mean(direct));
%
% COMPUTE 2D CASCADES
%
nt1=length(kj)/2;
direct=direct(1:nt1,1:nt1);
%
% seq1=-7.9e-4*1i*(f1_1+f2_1)*(f1_2+f2_2);
% seq2=-2.7e-3*1i*(conj(f1_1)+conj(f2_1))*(f1_2+f2_2);
% seq1=-7.9e-4*1i*(f1_1+f2_1+f3_1+f4_1)...
%     *(f1_2+f2_2+f3_2+f4_2);
% seq2=-2.7e-3*1i*(conj(f1_1)+conj(f2_1)+conj(f3_1)+conj(f4_1))...
%     *(f1_2+f2_2+f3_2+f4_2);

% par1=-8e-4*1i*(f1_2+f2_2)*(f1_p+f2_p);
% par2=-0.09*1i*(f1_2+f2_2)*(f1_p+f2_p);
% for j=1:nt1
% for jj=1:nt1
% % par1(jj,j)=-8e-4*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
% % par2(jj,j)=-0.09*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
% seq1(jj,j)=-1i*(f1(j)+f2(j)+f3(j)+f4(j))*(f1(jj)+f2(jj)+f3(jj)+f4(jj));%7.9e-4*
% seq2(jj,j)=-1i*(conj(f1(j))+conj(f2(j))+conj(f3(j))+conj(f4(j)))*(f1(jj)+f2(jj)+f3(jj)+f4(jj));%2.7e-3*
% par1(jj,j)=-1i*(f1(jj)+f2(jj)+f3(jj)+f4(jj))...
%     *(f1(j+jj)+f2(j+jj)+f3(j+jj)+f4(j+jj));%8e-4*
% par2(jj,j)=-1i*(f1(jj)+f2(jj)+f3(jj)+f4(jj))...
%     *(f1(j+jj)+f2(j+jj)+f3(j+jj)+f4(j+jj));%0.09*
% end
% end


for j=1:nt1
for jj=1:nt1
% par1(jj,j)=-8e-4*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
% par2(jj,j)=-0.09*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
seq1(jj,j)=-1i*(f1(j)+f2(j))*(f1(jj)+f2(jj));%7.9e-4*
seq2(jj,j)=-1i*(conj(f1(j))+conj(f2(j)))*(f1(jj)+f2(jj));%2.7e-3*
par1(jj,j)=-1i*(f1(jj)+f2(jj))...
    *(f1(j+jj)+f2(j+jj));%8e-4*
par2(jj,j)=-1i*(f1(jj)+f2(jj))...
    *(f1(j+jj)+f2(j+jj));%0.09*
end
end



seq=seq1(1:nt1,1:nt1)+seq2(1:nt1,1:nt1);
par=par1+par2;
cascade=seq+par;

% time=tau_1(1:nt1);
% wt=exp(-2*pi*1i*wvib*c*time);
% [wt1,wt2]=meshgrid(wt,wt);
% cas=abs(sum(sum(cascade.*wt1.*wt2)));
% dir=abs(sum(sum(direct.*wt1.*wt2)));
% ratio=cas/dir;
dw=1/nt1;
w=(-1/2:dw:(1/2-dw))/dt/c;
[~,iomega]=min(abs(w-wvib));

cas2=abs(fftshift(fft2(cascade)));
dir2=abs(fftshift(fft2(direct)));
cas=cas2(iomega,iomega);
dir=dir2(iomega,iomega);
ratio=cas/dir;
% cascade=cascade-mean(mean(cascade));
