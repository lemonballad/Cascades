%% This script calculates third order cascades for FSRS
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all
% Plot spectra from response functions
plot_flag=true;
% Speed of light cm/fs
c=2.998E-5;

%% Default material parameters
% Concentratin in mol L^-1
C=2E-4; % just a number
% unitless displacement
disp=0.35;%[0.35 0.35 1];
% electronic dephasing in cm^-1
gamma_eg=3000;%2010;
% vibrational dephasing in cm^-1
gamma_vib=20;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=6.2; % this is not correct, just a number
% refractive index
n_w_t=1.39; % also just a number
% Electronic energy gap origin
weg=22720;
% Vibrational frequencies
wvib=1600;
wvibs=[861 1599];
% Signal frequency
w_t=weg;
% Myoglobin parameters
if false
    [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters;
    wvib=wvib(1);disp=disp(1);
end
% PNA parameters
if true
     [~,gamma_eg,kappa,lambda,weg,translen,n_w_t,...
    wvibs,disp]=PNA_parameters('acetonitrile');
    wvibs=wvibs([1 5]);disp_pna=disp([1 5]);dts=1./wvibs/4/c;
end
% Set number of vibrational modes
nmode=1;%length(wvib);
% Set number of vibrational quanta to apportion
nquanta=5;
% Compute E(3):E(5) prefactor
prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t);

for iv=1:1
    %% Default laser paramters
    wvib=wvibs(iv);
    dt=dts(iv);
    w_L=weg;
    nt=100;
    %
    dws=-000:100:2000;nw=length(dws);
    disps=[ 0:0.1:0];nd=length(disps);
   for iw=1:nw
        for id=1:nd
            w_L=dws(iw)+weg-wvib;
            disp=disps(id);
            
            % Compute E(3):E(5) prefactor
            prefactor(id,iw,iv)=prefactor_3_5(l,C,mu_eg,n_w_t,w_L+wvib);
            
            %% Consolidate parameters
            % Consolidate material parameters into structure variable
            % parameters_material
            parameters_material.disp=disp;
            parameters_material.gamma_eg=gamma_eg;
            parameters_material.gamma_vib=gamma_vib;
            parameters_material.weg=weg;
            parameters_material.wvib=wvib;
            
            %     % Consolidate laser parameters into structure variable parameters_laser
            parameters_laser.dt=dt;
            parameters_laser.nt=nt;
            parameters_laser.w_L=w_L;
            
            %% Compute matrix of basis states, vibrational energies, and overlap integrals
            % Basis states and vibrational energies
            [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
            % Overlap integrals
            [ovlp_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);
            tic
            [rat,cas,dir] = cascadesub4_TCBETA(wviball_TC,nquanta,ovlp_TC,...
                parameters_laser,parameters_material,wvib);
            ratio(id,iw,iv)=rat*prefactor(id,iw,iv)/3e10;rcas(id,iw,iv)=cas;rdir(id,iw,iv)=dir;
            [toc iv iw id]
        end
    end
    % figure;semilogy(disps,ratio(:,(nw+1)/2));
    % figure;contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    % xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    % ylabel('Mode Displacement');
    % colormap jet;colorbar('Location','northoutside','fontsize',16)
    % set(gca,'linewidth',2,'fontsize',16);axis square;
    
    %% Spectra
    %% Default laser paramters
%     disp=disp_pna(iv);
%     wvib=wvibs(iv);
%     dt=8.25;
%     w_L=weg;
%     nt=75;
%     dw=1/nt;
%     w=(-1/2:dw:(1/2-dw))/dt/c;
%     [~,iomega]=min(abs(w-wvib));
    
%     %% Consolidate parameters
%     % Consolidate material parameters into structure variable
%     % parameters_material
%     parameters_material.disp=disp;
%     parameters_material.gamma_eg=gamma_eg;
%     parameters_material.gamma_vib=gamma_vib;
%     parameters_material.weg=weg;
%     parameters_material.wvib=wvib;
%     
%     %     % Consolidate laser parameters into structure variable parameters_laser
%     parameters_laser.dt=dt;
%     parameters_laser.nt=nt;
%     parameters_laser.w_L=weg;
    
    %% Compute matrix of basis states, vibrational energies, and overlap integrals
    % Basis states and vibrational energies
%     [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
%     % Overlap integrals
%     [ovlp_TC] = fcinfo_TC(base_TC,disp,nmode,nquanta);
%     tic
%     [~,~,cascade2d(:,:,iv),direct2d(:,:,iv)] = cascadesub4_TC2(wviball_TC,nquanta,ovlp_TC,...
%         parameters_laser,parameters_material);
%     toc
    
%     figure;
%     subplot(1,3,1);contour(w,w,abs(cascade2d(:,:,iv)).^0.5,50);
%     xlabel('\omega_1/2\pic (cm^{-1})');
%     ylabel('\omega_2/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',16)
%     set(gca,'linewidth',2,'fontsize',16);axis square;title('CASCADES')
%     subplot(1,3,2);contour(w,w,abs(direct2d(:,:,iv)).^0.5,50);
%     xlabel('\omega_1/2\pic (cm^{-1})');
%     ylabel('\omega_2/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',16)
%     set(gca,'linewidth',2,'fontsize',16);axis square;title('Direct')
%     subplot(1,3,3);contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
%     xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
%     ylabel('Mode Displacement');
%     colormap jet;colorbar('Location','northoutside','fontsize',16)
%     set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
    
end

if false%plot_flag
    seq=seq-seq(1,end);
    par=par-par(1,end);
    cascade2d=(cascade2d)-(cascade2d(1,100));
    direct2d=(direct2d)-(direct2d(1,100));
    figure;
    subplot(2,3,1);contour(w,w,abs(seq).^.5,25);colorbar;title('SEQUENTIAL CASCADES')
    subplot(2,3,2);contour(w,w,abs(par).^.5,25);colorbar;title('PARALLEL CASCADES')
    subplot(2,3,3);contour(w,w,abs(cascade2d).^.5,25);colorbar;title('CASCADES')
    subplot(2,3,4);contour(w,w,abs(direct2d).^.5,25);colorbar;title('DIRECT')
    subplot(2,3,6);contour(w,w,abs(cascade2d+direct2d).^.5,25);colorbar;title('DIRECT+CASCADES')
end


