%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
c=2.998E-5;%

%% Default material parameters
% Concentratin in mol L^-1
C=2E-4; % just a number
% unitless displacement
disp=0.35;%[0.35 0.35 1];
% electronic dephasing in cm^-1
gamma_eg=750;%2010;
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=8.8;%6.2; % this is not correct, just a number
% refractive index
n_w_t=1.39; % also just a number
% Electronic energy gap origin
weg=23250;%38000;
% Vibrational frequencies
wvibs=[670 1370];%[1370 1100 1500];
% Signal frequency
w_t=weg;
% Myoglobin parameters
if false
    [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters;
end
% PNA parameters
if false
     [~,gamma_eg,kappa,lambda,weg,translen,n_w_t,...
    wvib,disp]=PNA_parameters('acetonitrile');
    wvibs=wvib([1 5]);disp_pna=disp([1 5]);
end
% Set number of vibrational modes
nmode=1;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=7;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);

%% Default laser paramters
% Spectral width actinic pump
LAMBDA_ap=210;
% Spectral width raman pump
LAMBDA_rp=40;
% Frequency actinic pump
w_ap=weg;
% Frequency raman pump
w_rp=w_ap;
%
dt=100;%100;8.25;
%
w_aps=-000:50:1000;
disps=0:0.2:1.8;
nt=100;
nd=length(disps);
nw=length(w_aps);
% t=0:dt:nt*dt-dt;
% dw=1/nt;
% w=(-1/2:dw:(1/2-dw))/dt/c;w=w-w(1);
cascade_3=complex(zeros(nd,nw,2,nt));
direct=complex(zeros(nd,nw,2,nt));
rTC=zeros(nd,nw,2,'double');
nl=nw*nd*2;
for iv=1:2
    wvib=wvibs(iv);
    disp=disps(iv);
    for iw=1:nw
        for id=1:nd
            disp=disps(id);
            prefactor(iv)=prefactor_3_5(l,C,mu_eg,n_w_t,w_rp-wvib);
            w_ap=w_aps(iw)+weg-wvib/2;
            w_rp=w_ap;
%             w_t=w_rp-w;
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
            % Overlap integrals
            [fcall_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);
            
            %% Consolidate parameters
            % Consolidate material parameters into structure variable
            % parameters_material
            parameters_material.disp=disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
            
            % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.LAMBDA_ap=LAMBDA_ap;
            parameters_laser.LAMBDA_rp=LAMBDA_rp;
            parameters_laser.w_ap=w_ap;
            parameters_laser.w_rp=w_rp;
%             parameters_laser.w_t=w_t;
            parameters_laser.dt=dt;
            parameters_laser.nt=nt;
            
            
            %% Compute third order cascade response functions
            %         [d,s,e,~,~] = response2_Tuner(wviball_TC,fcall_TC,wvib,w_ap);
            %         rAM(id,iw)=3e13*abs(s(100))./abs(d(100))/3e10/2;
            tic
            [direct(id,iw,iv,:),cascade_3(id,iw,iv,:),~,~,~] = FSRS_TC_TuningNew(wviball_TC,fcall_TC,...
                3,parameters_material,parameters_laser);
            nl=nl-1;
            [toc nl]
                     
             rTC(id,iw,iv)=prefactor(iv)*abs(cascade_3(id,iw,iv,1))./abs(direct(id,iw,iv,1))/3e10;
        end
    end
end
