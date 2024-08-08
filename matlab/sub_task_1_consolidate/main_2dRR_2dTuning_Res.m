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
weg=22720;
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
% Time step for specific vibrational mode.
Res_2dRR_2d_dts=10;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;

% Initialize arrays to 0
Res_2dRR_2d_dws=-3000:50:3000;nw=length(Res_2dRR_2d_dws);
Res_2dRR_2d_dvibs=(50:50:3000);nd=length(Res_2dRR_2d_dvibs);
Res_2dRR_2d_nt=128;
Res_2dRR_2d_Cascade=complex(zeros(Res_2dRR_2d_nt,Res_2dRR_2d_nt,'double'));
Res_2dRR_2d_Direct=complex(zeros(Res_2dRR_2d_nt,Res_2dRR_2d_nt,'double'));
Res_2dRR_2d_Ratio=complex(zeros(nd,nw,'double'));
    %% Default laser paramters
    Disp=disp_pna;
    Res_2dRR_2d_dts=10;dt=Res_2dRR_2d_dts;
    matlabpool open 10
   for iw=1:nw
       tic
        for id=1:nd
            wvib=Res_2dRR_2d_dvibs(id);
            w_L=Res_2dRR_2d_dws(iw)+weg;
            % Compute E(3):E(5) prefactor
            prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_L);
            
            %% Consolidate parameters
            % Consolidate material parameters into structure variable 
            % parameters_material
            parameters_material.disp=Disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
           
            %     % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.dt=dt;
            parameters_laser.nt=Res_2dRR_2d_nt;
            parameters_laser.w_L=w_L;

            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);

            [rat,cas,dir] = cascade_2dRR_Res(wviball_TC,nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            Res_2dRR_2d_Ratio(id,iw)=rat*prefactor/3e10;
            Res_2dRR_2d_Cas(id,iw)=cas*prefactor/3e10;
            Res_2dRR_2d_Dir(id,iw)=dir;

%             if Res_2dRR_2d_dws(iw)==0&&Res_2dRR_2d_dvibs(id)==500
%                 if nquanta>3;nq=3;else nq=nquanta;end
%                 parameters_material.disp=disp_pna;
%                 [cas2,dir2] = cascade_2dRR_ResSpec(wviball_TC,nq,ovlp_TC,...
%                     parameters_laser,parameters_material);
%                 Res_2dRR_2d_Cascade(:,:)=prefactor/3e10*cas2;
%                 Res_2dRR_2d_Direct(:,:)=dir2;
%             end
        end
        disp([num2str(toc) ' s, dw = ' num2str(Res_2dRR_2d_dws(iw))]);
   end
    save('2dRR2d','Res_2dRR_2d_Ratio','Res_2dRR_2d_Cascade',...
        'Res_2dRR_2d_Direct','Res_2dRR_2d_dts','Res_2dRR_2d_nt',...
        'Res_2dRR_2d_dws','Res_2dRR_2d_dvibs','Res_2dRR_2d_Cas','Res_2dRR_2d_Dir');
matlabpool close
    'done'

