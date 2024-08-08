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
C=2E-4; % 
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length in m
l=2.2E-4; % 

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
FSRS_Res_dts=10;1./wvibs/2/c;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;

%% Default laser paramters
% Spectral width actinic pump
LAMBDA_ap=210;
% Spectral width raman pump
LAMBDA_rp=40;
% Number of time steps
FSRS_Res_nt=128;

% Range of laser detuning
FSRS_Res_dws=-1000:250:3000;
% Length of laser detuning range
FSRS_Res_nw=length(FSRS_Res_dws);
% Range of mode displacement
FSRS_Res_disps=0.2:0.2:3;
% Length of mode displacement range
nd=length(FSRS_Res_disps);
% Initialize storage arrays to correct size of 0s
FSRS_Res_Cascade=complex(zeros(FSRS_Res_nt,FSRS_Res_nt,2,'double'));
FSRS_Res_Direct=complex(zeros(FSRS_Res_nt,FSRS_Res_nt,2,'double'));
FSRS_Res_Ratio=complex(zeros(nd,FSRS_Res_nw,2,'double'));
matlabpool open 10
% Loop over vibrational modes
for iv=1:2
    % Set vibrational mode for specific loop iteration
    wvib=wvibs(iv);
    % Set time step for specific loop iteration
    dt=FSRS_Res_dts;
    
    % Loop of laser detuning
    for iw=1:length(FSRS_Res_dws)
        % Set center frequency of actinic pump
        w_ap=FSRS_Res_dws(iw)+weg;
        % Set center frequency of Raman pump
        w_rp=w_ap;
        % Calculate prefactor
        prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_rp-wvib);
        tic
        
        % Loop over mode displacements
        for id=1:length(FSRS_Res_disps)
            % Set mode displacement for specific loop iteration
            Disp=disp_pna*FSRS_Res_disps(id);
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
            % Overlap integrals
           [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);
            
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
            parameters_laser.nt=FSRS_Res_nt;
            
            %% Compute third order cascade response functions
%             tic
            % Get rat=Ratio, cas=Cascade, and dir=Direct for all resonant
            % FSRS
            [rat,cas,dir] = cascade_FSRS_Res(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            FSRS_Res_Ratio(id,iw,iv)=prefactor/3e10*rat;
            rcas(id,iw,iv)=prefactor/3e10*cas;rdir(id,iw,iv)=dir;

            % Get full spectra for 0 detuning and mode displacements for
            % specific vibrational modes
            if (FSRS_Res_dws(iw)==0)&&FSRS_Res_disps(id)==0.1
                % Keep number of quanta for spectra at or below 5 to reduce
                % computational cost
                if nquanta>=3,nq=3;else nq=nquanta;end
                %% Compute matrix of basis states, vibrational energies, and overlap integrals
                % Basis states and vibrational energies
                [base_TC,wviball_TC]=basis_TC(nmode,nq,wvibs);
                % Overlap integrals
                [ovlp_TC] = fcinfo_TC(base_TC,disp_pna,nmode,nq);

                % Define mode displacement
                parameters_material.disp=disp_pna(iv);
                % Calculate spectra for cascade and direct signal
                [FSRS_Res_Cascade(:,:,iv),FSRS_Res_Direct(:,:,iv)] = cascade_FSRS_ResSpec(...
                    wviball_TC,2*nq,ovlp_TC,parameters_laser,parameters_material);
                FSRS_Res_Cascade(:,:,iv)=prefactor/3e10*FSRS_Res_Cascade(:,:,iv);
            end
            
        end
        disp([num2str(toc) ' s, dw = ' num2str(FSRS_Res_dws(iw))]);
    end
end
% Save data for plotting
save('FSRS','FSRS_Res_Ratio','FSRS_Res_Cascade','FSRS_Res_Direct',...
    'FSRS_Res_dts','FSRS_Res_nt','FSRS_Res_dws','FSRS_Res_disps','rcas','rdir');
matlabpool close
'done'