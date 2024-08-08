%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
% Define constants
c=2.998E-5;% Speed of light cm/fs
q_e=1.60218e-19; % Fundamental charge of an electron in C
Dconv=3.33564e-30; % Conversion factor to D
kT=200; % kT in cm^-1

%% Default material parameters
% Concentratin in mol L^-1
C=2E-4; % just a number
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% Electronic energy gap origin
weg=20000;%38000;
% Signal frequency
w_t=weg;

% PNA parameters
[~,~,kappa,lambda,weg,translen,n_w_t,...
    wvib,Disp]=PNA_parameters('acetonitrile');
% Vibrational modes to be examined in cm^-1
wvibs=wvib([1 3]);
% Unitless mode displacements
disp_pna=Disp([1 3]);
% Electronic dephasing in cm^-1
gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
% Transition dipole in D
mu_eg=translen*1e-10*q_e/Dconv;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;
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
FSRS_Res_2d_nt=128;

FSRS_Res_2d_dws=-3000:50:3000;
FSRS_Res_2d_dvibs=50:50:3000;
nd=length(FSRS_Res_2d_dvibs);
FSRS_Res_2d_nw=length(FSRS_Res_2d_dws);
FSRS_Res_2d_Cascade=complex(zeros(FSRS_Res_2d_nt,FSRS_Res_2d_nt,'double'));
FSRS_Res_2d_Direct=complex(zeros(FSRS_Res_2d_nt,FSRS_Res_2d_nt,'double'));
FSRS_Res_2d_Ratio=complex(zeros(nd,FSRS_Res_2d_nw,'double'));
FSRS_Res_2d_dts=10;
dt=FSRS_Res_2d_dts;
Disp=disp_pna;
matlabpool open 10
    for iw=1:length(FSRS_Res_2d_dws)
        tic
        for id=1:length(FSRS_Res_2d_dvibs)
            wvib=FSRS_Res_2d_dvibs(id);
            w_ap=FSRS_Res_2d_dws(iw)+weg;
            w_rp=w_ap;
            prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_rp-wvib);
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,disp_pna,nmode,nquanta);
            
            %% Consolidate parameters
            % Consolidate material parameters into structure variable
            % parameters_material
            parameters_material.disp=Disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
            
            % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.LAMBDA_ap=LAMBDA_ap;
            parameters_laser.LAMBDA_rp=LAMBDA_rp;
            parameters_laser.w_ap=w_ap;
            parameters_laser.w_rp=w_rp;
            parameters_laser.dt=dt;
            parameters_laser.nt=FSRS_Res_2d_nt;
            
            %% Compute third order cascade response functions
            [rat,cas,dir] = cascade_FSRS_Res(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            FSRS_Res_2d_Ratio(id,iw)=prefactor/3e10*rat;
            FSRS_Res_2d_Cas(id,iw)=prefactor/3e10*cas;
            FSRS_Res_2d_Dir(id,iw)=dir;
            
%             if FSRS_Res_2d_dws(iw)==0&&wvib==200
%                 if nquanta>=3,nq=3;else nq=nquanta;end
%             %% Compute matrix of basis states, vibrational energies, and overlap integrals
%             % Basis states and vibrational energies
%             [base_TC,wviball_TC]=basis_TC(nmode,nq,wvibs);
%             % Overlap integrals
%             [ovlp_TC] = fcinfo_TC(base_TC,disp_pna,nmode,nq);
%                 [FSRS_Res_2d_Cascade(:,:),FSRS_Res_2d_Direct(:,:)] = cascade_FSRS_ResSpec(...
%                     wviball_TC,2*nq,ovlp_TC,parameters_laser,parameters_material);
%                 FSRS_Res_2d_Cascade(:,:)=prefactor/3e10*FSRS_Res_2d_Cascade(:,:);
%             end
            
        end
        disp([num2str(toc) ' s, dw = ' num2str(FSRS_Res_2d_dws(iw))]);
    end

save('FSRS2d','FSRS_Res_2d_Ratio','FSRS_Res_2d_Cascade','FSRS_Res_2d_Direct',...
    'FSRS_Res_2d_dts','FSRS_Res_2d_nt','FSRS_Res_2d_dws','FSRS_Res_2d_dvibs',...
    'FSRS_Res_2d_Cas','FSRS_Res_2d_Dir');
matlabpool close
'done'