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
Res_2dRR_dts=5;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=3;

% Initialize arrays to 0
Res_2dRR_dws=-1000:250:3000;nw=length(Res_2dRR_dws);
Res_2dRR_disps=(0.2:0.2:3);nd=length(Res_2dRR_disps);
Res_2dRR_nt=64;
Res_2dRR_Cascade=complex(zeros(Res_2dRR_nt,Res_2dRR_nt,2,'double'));
Res_2dRR_Direct=complex(zeros(Res_2dRR_nt,Res_2dRR_nt,2,'double'));
Res_2dRR_Ratio=complex(zeros(nd,nw,2,'double'));
% matlabpool open 10
for iv=1:2
    %% Default laser paramters
    wvib=wvibs(iv);
    dt=Res_2dRR_dts;
   for iw=1:nw
       tic
        for id=1:nd
            w_L=Res_2dRR_dws(iw)+weg;
            Disp=disp_pna*Res_2dRR_disps(id);
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
            parameters_laser.nt=Res_2dRR_nt;
            parameters_laser.w_L=w_L;

            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);

            [rat,cas,dir] = cascade_2dRR_Res(wviball_TC,nmode*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            Res_2dRR_Ratio(id,iw,iv)=rat*prefactor/3e10;
            Res_2dRR_Cas(id,iw,iv)=cas*prefactor/3e10;
            Res_2dRR_Dir(id,iw,iv)=dir;

            if Res_2dRR_dws(iw)==0&&Res_2dRR_disps(id)==0.2
                if nquanta>3,nq=3;else nq=nquanta;end
                %% Compute matrix of basis states, vibrational energies, and overlap integrals
                % Basis states and vibrational energies
                [base_TC,wviball_TC]=basis_TC(nmode,nq,wvibs);
                % Overlap integrals
                [ovlp_TC] = fcinfo_TC(base_TC,disp_pna,nmode,nq);
                parameters_material.disp=disp_pna(iv);
                [cas2,dir2] = cascade_2dRR_ResSpec(wviball_TC,nmode*nq,ovlp_TC,...
                    parameters_laser,parameters_material);
                Res_2dRR_Cascade(:,:,iv)=prefactor/3e10*cas2;
                Res_2dRR_Direct(:,:,iv)=dir2;
            end
        end
        disp([num2str(toc) ' s, dw = ' num2str(Res_2dRR_dws(iw))]);
   end
end
    save('2dRRtest','Res_2dRR_Ratio','Res_2dRR_Cascade',...
        'Res_2dRR_Direct','Res_2dRR_dts','Res_2dRR_nt',...
        'Res_2dRR_dws','Res_2dRR_disps','Res_2dRR_Cas','Res_2dRR_Dir');
% matlabpool close
'done'
