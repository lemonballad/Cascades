function [ratio,cas,dir] = cascadesub4_TC3(E_vib,nquanta,ovlp,...
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
kj=1:nt;
taud=(kj-1)*dt;
kjf=1:nt*2;
tau=(kjf-1)*dt;
nt1=nt*2;

% Initialize terms of response function to 0
r1=complex(zeros(nt,nt,iq,'double'));
r2=r1;r3=r1;r4=r1;
r5=r1;r6=r1;r7=r1;r8=r1;
r9=r1;r10=r1;r11=r1;r12=r1;
r13=r1;r14=r1;r15=r1;r16=r1;

f1=complex(zeros(1,nt1,iq,'double'));
f2=f1;f3=f1;f4=f1;


[tau_1,tau_2]=meshgrid(taud);
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
                if k~=m
                    f1(:,:,n)=f1(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,m)...
                        *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                        .*exp((-1i*w(k,m)*r2w-gamma_vib*c)*tau)...
                        *1/(w_L-weg-w(l,m)+1i*gamma_eg);
                    %
                    f2(:,:,n)=f2(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,k)...
                        *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                        .*exp((-1i*w(m,k)*r2w-gamma_vib*c)*tau)...
                        *1/(w_L-weg-w(l,k)+1i*gamma_eg);
                end
                %
                if k~=n
                    f3(:,:,n)=f3(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,l)...
                        *1./(w_L-weg-w(n,m)+1i*gamma_eg)...
                        .*exp((-1i*w(n,k)*r2w-gamma_vib*c)*tau)...
                        *1./(w_L-weg-w(n,l)+1i*gamma_eg);
                    %
                    f4(:,:,n)=f4(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,l)...
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
                        if m~=u && m~=k
                            % Auxillary function #1 for fifth order signal
                            r1(:,:,n)=r1(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,u)*ovlp(v,m)...
                                *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                *1./(w_L-weg-w(l,m)+1i*gamma_eg)...
                                .*exp((-1i*w(u,m)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,m)+1i*gamma_eg);
                            
                        end
                        %
                        if k~=u && m~=k
                            % Auxillary function #2 for fifth order signal
                            r2(:,:,n)=r2(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,k)*ovlp(v,u)...
                                *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                *1./(-w_L+weg-w(k,l)+1i*gamma_eg)...
                                .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                        end
                        %
                        if m~=u && m~=k
                            % Auxillary function #3 for fifth order signal
                            r3(:,:,n)=r3(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(l,u)*ovlp(v,m)*ovlp(v,u)...
                                *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                *1./(-w_L+weg-w(m,l)+1i*gamma_eg)...
                                .*exp((-1i*w(m,u)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                        end
                        %
                        if k~=u && m~=k
                            % Auxillary function #4 for fifth order signal
                            r4(:,:,n)=r4(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(l,u)*ovlp(v,u)*ovlp(v,k)...
                                *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                *1./(w_L-weg-w(l,k)+1i*gamma_eg)...
                                .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
                                *1./(w_L-weg-w(v,k)+1i*gamma_eg);
                        end
                        %
                        %
                        if true
                            if k~=u && n~=k
                                % Auxillary function #5 for fifth order signal
                                r5(:,:,n)=r5(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(k,v)*ovlp(u,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,k)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,k)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                            end
                            %
                            if n~=u && n~=k
                                % Auxillary function #6 for fifth order signal
                                r6(:,:,n)=r6(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(u,v)*ovlp(n,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(n,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(n,v)+1i*gamma_eg);
                                %
                                % Auxillary function #7 for fifth order signal
                                r7(:,:,n)=r7(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(u,l)*ovlp(n,v)*ovlp(u,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,n)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                            end
                            %
                            if k~=u && n~=k
                                % Auxillary function #8 for fifth order signal
                                r8(:,:,n)=r8(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(u,l)*ovlp(u,v)*ovlp(k,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(k,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(k,v)+1i*gamma_eg);
                            end
                            %
                            if l~=u && m~=k
                                % Auxillary function #9 for fifth order signal
                                r9(:,:,n)=r9(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(u,v)*ovlp(l,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(l,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(l,v)+1i*gamma_eg);
                                %
                                % Auxillary function #10 for fifth order signal
                                r10(:,:,n)=r10(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(l,v)*ovlp(u,v)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,m)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(k,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                                %
                                % Auxillary function #11 for fifth order signal
                                r11(:,:,n)=r11(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,m)*ovlp(u,k)*ovlp(u,v)*ovlp(l,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(l,k)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(l,v)+1i*gamma_eg);
                                %
                                % Auxillary function #12 for fifth order signal
                                r12(:,:,n)=r12(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(n,k)*ovlp(l,k)*ovlp(u,m)*ovlp(l,v)*ovlp(u,v)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(m,k)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(m,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(u,v)+1i*gamma_eg);
                            end
                            if l~=u && n~=k
                                %
                                % Auxillary function #13 for fifth order signal
                                r13(:,:,n)=r13(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,l)*ovlp(v,u)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,k)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                                %
                                % Auxillary function #14 for fifth order signal
                                r14(:,:,n)=r14(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,u)*ovlp(v,l)...
                                    *1/(w_L-weg-w(n,m)+1i*gamma_eg)...
                                    .*exp((-1i*w(n,k)*r2w-damp)*tau_1 )...
                                    *1./(w_L-weg-w(n,l)+1i*gamma_eg)...
                                    .*exp((-1i*w(u,l)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,l)+1i*gamma_eg);
                                %
                                % Auxillary function #15 for fifth order signal
                                r15(:,:,n)=r15(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(k,l)*ovlp(n,u)*ovlp(v,l)*ovlp(v,u)...
                                    *1/(-w_L+weg-w(m,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(k,n)*r2w-damp)*tau_1 )...
                                    *1./(-w_L+weg-w(l,n)+1i*gamma_eg)...
                                    .*exp((-1i*w(l,u)*r2w-gamma_vib*c)*tau_2)...
                                    *1./(w_L-weg-w(v,u)+1i*gamma_eg);
                                %
                                % Auxillary function #16 for fifth order signal
                                r16(:,:,n)=r16(:,:,n)+boltz_pop(m)*ovlp(n,m)*ovlp(k,m)*ovlp(n,l)*ovlp(k,u)*ovlp(v,u)*ovlp(v,l)...
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
    f1=sum(f1,3);f2=sum(f2,3);f3=sum(f3,3);f4=sum(f4,3);
    r1=sum(r1,3);r2=sum(r2,3);r3=sum(r3,3);r4=sum(r4,3);
    r5=sum(r5,3);r6=sum(r6,3);r7=sum(r7,3);r8=sum(r8,3);
    r9=sum(r9,3);r10=sum(r10,3);r11=sum(r11,3);r12=sum(r12,3);
    r13=sum(r13,3);r14=sum(r14,3);r15=sum(r15,3);r16=sum(r16,3);
end % End loop over initial state

%
% direct=-0.09*(r1+r2+r3+r4);
direct=-((r1+r2+r3+r4)+(r5+r6+r7+r8)...
    +(r9+r10+r11+r12)+(r13+r14+r15+r16));%0.09*
% direct=direct-mean(mean(direct));
%
% COMPUTE 2D CASCADES
%
% nt1=length(kjf)/2;
nt1=nt;

% direct=direct(1:nt1,1:nt1);
%
% seq1=-7.9e-4*1i*(f1_1+f2_1)*(f1_2+f2_2);
% seq2=-2.7e-3*1i*(conj(f1_1)+conj(f2_1))*(f1_2+f2_2);
% seq1=-7.9e-4*1i*(f1_1+f2_1+f3_1+f4_1)...
%     *(f1_2+f2_2+f3_2+f4_2);
% seq2=-2.7e-3*1i*(conj(f1_1)+conj(f2_1)+conj(f3_1)+conj(f4_1))...
%     *(f1_2+f2_2+f3_2+f4_2);

% par1=-8e-4*1i*(f1_2+f2_2)*(f1_p+f2_p);
% par2=-0.09*1i*(f1_2+f2_2)*(f1_p+f2_p);
for j=1:nt1
for jj=1:nt1
% par1(jj,j)=-8e-4*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
% par2(jj,j)=-0.09*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
seq1(jj,j)=-1i*(f1(j)+f2(j)+f3(j)+f4(j))*(f1(jj)+f2(jj)+f3(jj)+f4(jj));%7.9e-4*
seq2(jj,j)=-1i*(conj(f1(j))+conj(f2(j))+conj(f3(j))+conj(f4(j)))*(f1(jj)+f2(jj)+f3(jj)+f4(jj));%2.7e-3*
par1(jj,j)=-1i*(f1(jj)+f2(jj)+f3(jj)+f4(jj))...
    *(f1(j+jj)+f2(j+jj)+f3(j+jj)+f4(j+jj));%8e-4*
par2(jj,j)=-1i*(f1(jj)+f2(jj)+f3(jj)+f4(jj))...
    *(f1(j+jj)+f2(j+jj)+f3(j+jj)+f4(j+jj));%0.09*
end
end


% for j=1:nt1
% for jj=1:nt1
% % par1(jj,j)=-8e-4*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
% % par2(jj,j)=-0.09*1i*(f1(jj)+f2(jj))*(f1(j+jj)+f2(j+jj));
% seq1(jj,j)=-1i*(f1(j)+f2(j))*(f1(jj)+f2(jj));%7.9e-4*
% seq2(jj,j)=-1i*(conj(f1(j))+conj(f2(j)))*(f1(jj)+f2(jj));%2.7e-3*
% par1(jj,j)=-1i*(f1(jj)+f2(jj))...
%     *(f1(j+jj)+f2(j+jj));%8e-4*
% par2(jj,j)=-1i*(f1(jj)+f2(jj))...
%     *(f1(j+jj)+f2(j+jj));%0.09*
% end
% end



seq=seq1(1:nt1,1:nt1)+seq2(1:nt1,1:nt1);
par=par1+par2;
cascade=seq+par;

% time=tau_1(1:nt1);
% wt=exp(-2*pi*1i*omega*c*time);
% [wt1,wt2]=meshgrid(wt,wt);
dw=1/nt1;
w=(-1/2:dw:(1/2-dw))/dt/c;
[~,iomega]=min(abs(w-wvib));

cas2=abs(fftshift(fft2(cascade)));
dir2=abs(fftshift(fft2(direct)));
cas=cas2(iomega,iomega);
dir=dir2(iomega,iomega);
ratio=cas/dir;
% cascade=cascade-mean(mean(cascade));

%
% FOURIER TRANSFORM tau_1 & tau_2
%
% seq=fftshift(fft2(seq));
% par=fftshift(fft2(par));
% cascade2d=fftshift(fft2(cascade));
% direct2d=fftshift(fft2(direct));
% ratio=abs(cascade2d(omega,omega)/abs(direct2d(omega,omega);
% test=direct;


% figure;
% subplot(1,2,1);contour(w,w,abs(cas2),50);colorbar
% subplot(1,2,2);contour(w,w,abs(dir2),50);colorbar
