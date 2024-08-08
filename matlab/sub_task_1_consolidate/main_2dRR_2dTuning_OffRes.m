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
% Solvent Raman differential cross section prefactor in m^2
K=7.62e-31; % Methanol
% Solvent density in g/m^3
rho_Solv=792000; % Methanol
% Solvent molar mass in g/mol
mm_Solv=32.04; % Methanol
% Solvent vibrational mode in cm^-1
wsolv=1599; % Methanol
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% refractive index Solute
n_w_t_Sol=1.39; % also just a number
% refractive index Solvent
n_w_t_Solv=1.34; % Methanol
% Electronic energy gap origin Solute in cm^-1
weg=22720;
% Electronic energy gap origin Solvent in cm^-1
weg_Solv=153100; % Methanol
% Signal frequency
w_t=weg;
% PNA parameters
[~,~,kappa,lambda,weg,translen,n_w_t,...
    wvib,Disp]=PNA_parameters('methanol');
% Vibrational modes to be examined in cm^-1
wvibs=wvib([1 3]);
% Unitless mode displacements
disp_pna=Disp([1 3]);
% Electronic dephasing in cm^-1
gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
% Transition dipole in D
mu_eg=translen*1e-10*q_e/Dconv;
% Time step for specific vibrational mode.
OffRes_2dRR_2d_dts=10;

% Set number of vibrational modes
nmode=2;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;

w_L=weg;
% Compute Solvent polarizability squared in A^4 s^8 kg^-2
alpha2=polarizability(K,weg_Solv,w_L);
% Compute Solvent concentration in mol m^-3
C_Solv=rho_Solv/mm_Solv;
C_Sol=2e-4;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t);
prefactorres=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_t);
OffRes_2dRR_2d_nt=128;
OffRes_2dRR_2d_dws=-3000:50:3000;nw=length(OffRes_2dRR_2d_dws);
OffRes_2dRR_2d_dvibs=(500:100:3000);nd=length(OffRes_2dRR_2d_dvibs);
OffRes_2dRR_2d_Cascade=complex(zeros(OffRes_2dRR_2d_nt,OffRes_2dRR_2d_nt,2,'double'));
OffRes_2dRR_2d_Direct=complex(zeros(OffRes_2dRR_2d_nt,OffRes_2dRR_2d_nt,2,'double'));
OffRes_2dRR_2d_Ratio=complex(zeros(nd,nw,2,'double'));
matlabpool open 10
for iv=1:2
    %% Default laser paramters
    wvib=wvibs(iv);
    Disp=disp_pna;
    dt=OffRes_2dRR_2d_dts;
   for iw=1:nw
       tic
        for id=1:nd
            w_L=OffRes_2dRR_2d_dws(iw)+weg;
            wsolv=OffRes_2dRR_2d_dvibs(id);
            
            % Compute E(3):E(5) prefactor
            alpha2=polarizability(K,weg_Solv,w_L);
            prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t)/4/pi^2;
            prefactorres=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_t);

            %% Consolidate parameters
            % Consolidate material parameters into structure variable
            % parameters_material
            parameters_material.disp=Disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
            parameters_material.wsolv=wsolv;
            
            %     % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.dt=dt;
            parameters_laser.nt=OffRes_2dRR_2d_nt;
            parameters_laser.w_L=w_L;
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);
            [rat,cas,dir] = cascade_2dRR_OffRes(wviball_TC,2*nquanta,ovlp_TC,...
                parameters_laser,parameters_material);
            OffRes_2dRR_2d_Ratio(id,iw,iv)=rat*prefactor*3e10;
            OffRes_2dRR_2d_Cas(id,iw,iv)=cas*prefactor*3e10;
            OffRes_2dRR_2d_Dir(id,iw,iv)=dir;
%             if OffRes_2dRR_2d_dws(iw)==0&&OffRes_2dRR_2d_dvibs(id)==0
%                 if nquanta>=5,nq=5;else nq=nquanta;end
%                 parameters_material.disp=disp_pna(iv);
%                 [OffRes_2dRR_2d_Cascade(:,:,iv),OffRes_2dRR_2d_Direct(:,:,iv)] = cascade_2dRR_OffResSpec(wviball_TC,nq,ovlp_TC,...
%                     parameters_laser,parameters_material);
%                 OffRes_2dRR_2d_Cascade(:,:,iv)=prefactor(id,iw,iv)*3e10*OffRes_2dRR_2d_Cascade(:,:,iv);
%            end
        end
        disp([num2str(toc) ' s, dw = ' num2str(OffRes_2dRR_2d_dws(iw))]);
    end
end
save('2dOffRR2d','OffRes_2dRR_2d_Ratio','OffRes_2dRR_2d_Cascade',...
    'OffRes_2dRR_2d_Direct','OffRes_2dRR_2d_dts','OffRes_2dRR_2d_nt',...
    'OffRes_2dRR_2d_dws','OffRes_2dRR_2d_dvibs','OffRes_2dRR_2d_Cas','OffRes_2dRR_2d_Dir');
matlabpool close
'done'