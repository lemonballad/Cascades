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
e0=8.85418782e-12;% Permitivity of free space in A^2 s^4 kg^-1 m^-3
solvents={'cyclohexane','1,4-dioxane','dichloromethane','acetonitrile','methanol'};

%% Default material parameters
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
XSect_FSRS_nt=128;
XSect_FSRS_Ratio=complex(zeros(5,2,2,'double'));
for isolv=1:5
    switch isolv
        case 1
            solvent='cyclohexane';
        case 2
            solvent='1,4-dioxane';
        case 3
            solvent='dichloromethane';
        case 4
            solvent='acetonitrile';
        case 5
            solvent='methanol';
        otherwise
            error('Error: This is not a valid solvent');
    end
    solvent
    % PNA parameters
    [~,~,kappa,lambda,weg,translen,n_w_t_Sol,...
        wvib,Disp,K,rho_Solv,mm_Solv,wsolv,n_w_t_Solv,...
        weg_Solv]=DEFENSE_PNA_parameters(solvent);
    % Signal frequency
    w_t=weg;
    % Vibrational modes to be examined in cm^-1
    wvibs=wvib([1 3]);
    % Unitless mode displacements
    disp_pna=Disp([1 3]);
    % Electronic dephasing in cm^-1
    gamma_eg=2*sqrt(2*log(2))*sqrt(2*kT*lambda);
    % Transition dipole in D
    mu_eg=translen*1e-10*q_e/Dconv;
    % Time step for specific vibrational mode.
    XSect_FSRS_dts=10;
    
    % Set number of vibrational modes
    nmode=2;%length(wvib);
    % Set number of vibrational quanta to apportion
    nquanta=5;
    
    % Compute Solvent concentration in mol m^-3
    C_Solv=rho_Solv/mm_Solv;
    C_Sol=2e-4;
    dt=XSect_FSRS_dts;
    w_L=weg;
    Disp=disp_pna;
    %% Default laser paramters
    % Spectral width actinic pump
    LAMBDA_ap=210;
    % Spectral width raman pump
    LAMBDA_rp=40;
    % Frequency actinic pump
    w_ap=weg;
    % Frequency raman pump
    w_rp=w_ap;

    % Compute E(3):E(5) prefactor
    alpha2=polarizability(K,weg_Solv,w_ap);
    prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_ap)/4/pi^2;
    prefactorres=prefactor_3_5(l,C_Sol,mu_eg,n_w_t_Sol,w_ap);
    Ks(isolv)=alpha2;%/e0^2*w_ap^4*10e8*10e20;
    Wsolvs(isolv)=wsolv;
    Rhos(isolv)=rho_Solv;
    for iv=1:2
        %% Default laser paramters
        wvib=wvibs(iv);
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
        parameters_laser.nt=XSect_FSRS_nt;
        parameters_laser.LAMBDA_ap=LAMBDA_ap;
        parameters_laser.LAMBDA_rp=LAMBDA_rp;
        parameters_laser.w_ap=w_ap;
        parameters_laser.w_rp=w_rp;

        %% Compute matrix of basis states, vibrational energies, and overlap integrals
        % Basis states and vibrational energies
        [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvibs);
        % Overlap integrals
        [ovlp_TC] = fcinfo_TC(base_TC,Disp,nmode,nquanta);
        
        [ratO,~,~] = cascade_FSRS_OffRes(wviball_TC,2*nquanta,ovlp_TC,...
            parameters_laser,parameters_material);
        XSect_FSRS_Ratio(isolv,iv,1)=ratO*prefactor*3e10;
        [ratR,~,~] = cascade_FSRS_Res(wviball_TC,2*nquanta,ovlp_TC,...
            parameters_laser,parameters_material);
        XSect_FSRS_Ratio(isolv,iv,2)=ratR*prefactorres/3e10;
        
        
    end
end
save('XSectFSRS','XSect_FSRS_Ratio');
'done'

[Ks,iKs]=sort(Ks);
figure;plot(Ks,0*XSect_FSRS_Ratio(iKs,2,2)...
    +XSect_FSRS_Ratio(iKs,2,1),'b^','MarkerFaceColor',[0 0 1]);
ylabel('RATIO');hold on;%axis square
plot(Ks,0*XSect_FSRS_Ratio(iKs,1,2)...
    +XSect_FSRS_Ratio(iKs,1,1),'r^','MarkerFaceColor',[1 0 0])
set(gca,'XTick',Ks,'XTickLabel',solvents(iKs),...
    'XTickLabelRotation',45,'FontSize',14,'FontWeight','Bold');
legend('~ 1330 cm{-1}','~ 860 cm{-1}');
hold off

% [Wsolvs,iWsolvs]=sort(Wsolvs);
% figure; plot(Wsolvs,0*XSect_FSRS_Ratio(iWsolvs,1,2)...
%     +XSect_FSRS_Ratio(iWsolvs,1,1),'r^','MarkerFaceColor',[1 0 0])
% hold on
% plot(Wsolvs,0*XSect_FSRS_Ratio(iWsolvs,2,2)...
%     +XSect_FSRS_Ratio(iWsolvs,2,1),'b^','MarkerFaceColor',[0 0 1])
% set(gca,'XTick',Wsolvs,'XTickLabel',solvents(iWsolvs),...
%     'XTickLabelRotation',45,'FontSize',14,'FontWeight','Bold');
% hold off
% 
% [Rhos,iRhos]=sort(Rhos);
% figure; plot(Rhos,0*XSect_FSRS_Ratio(iRhos,1,2)...
%     +XSect_FSRS_Ratio(iRhos,1,1),'r^','MarkerFaceColor',[1 0 0])
% hold on
% plot(Rhos,0*XSect_FSRS_Ratio(iRhos,2,2)...
%     +XSect_FSRS_Ratio(iRhos,2,1),'b^','MarkerFaceColor',[0 0 1])
% set(gca,'XTick',Rhos,'XTickLabel',solvents(iRhos),...
%     'XTickLabelRotation',45,'FontSize',14,'FontWeight','Bold');
% hold off