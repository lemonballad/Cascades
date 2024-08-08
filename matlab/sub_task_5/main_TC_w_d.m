%% This script calculates third order cascades
% Clear workspace
clear
% Clear command window
clc
% Clear open figures
% close all

%% Default material parameters
% Concentratin in mol L^-1
C=2E-4; % just a number
% unitless displacement
disp=0.35;%[0.35 0.35 1];
% electronic dephasing in cm^-1
gamma_eg=1000;
% vibrational dephasing in cm^-1
gamma_vib=10;
% path length
l=2.2E-4; % this seems right for the path length
% Transition dipole
mu_eg=8.8; % this is not correct, just a number
% refractive index
n_w_t=1.39; % also just a number
% Electronic energy gap origin
weg=20000;
% Vibrational frequencies
wvib=1370;%[1370 1100 1500];
% Signal frequency
w_t=weg;
% Myoglobin parameters
if false
    [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters;
end
% PNA parameters
if false
     [inhomog,gamma_0,kappa,lambda,weg,translen,refrac_index,...
    wvib,disp]=PNA_parameters('methanol');
end
% Set number of vibrational modes
nmode=length(wvib);
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

f_c=(sinc(-18*.022/2)*sinc(-197*.022/2)*exp(1i*(-18-197)*.022/2)...
    +sinc(-412*.022/2)*sinc(196*.022/2)*exp(1i*(-412+196)*.022/2));
f_d=(sinc(-215*.022/2)*exp(1i*-215*.022/2));

wvibs=500:100:1700;
disps=0:0.1:1;
for iw=1:length(wvibs)
    for id=1:length(disps)
        wvib=wvibs(iw);
        disp=disps(id);
        prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_rp+wvib);
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
        
        %% Compute matrix of basis states, vibrational energies, and overlap integrals
        % Basis states and vibrational energies
        [base_TC,wviball_TC]=basis_TC(nmode,nquanta,wvib);
        % Overlap integrals
        [fcall_tc] = fcinfo_TC(base_TC,disp,nmode,nquanta);
        
        %% Compute third order cascade response functions
        % [direct,seq,e3,tau,wrs] = response2_TC(wviball_TC,fcall_tc,parameters);
        [direct,cascade_3,e3,tau,wrs] = FSRS_TC_w_d(wviball_TC,fcall_tc,...
            4,parameters_material,parameters_laser);
        
%         ratio(iw,id)=3e13*abs(trapz(cascade_3))/abs(trapz(direct))/3e10/2;
        ratio(iw,id)=3e13*abs(f_c*cascade_3(1))./abs(f_d*direct(1))/3e10/2;
    end
end



%% Plot
if false
     ratio(isnan(ratio))=0;
    figure;contour(disps,wvibs,abs(ratio),35);%     axis square;
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);
    caxis([0.0 0.1]);
    
end

if false
    figure;contourf(wrs,tau,real(direct),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Real Component';'Direct Fifth Order'});
    
    figure;contourf(wrs,tau,imag(direct),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Imagninary Component';'Direct Fifth Order'});
    
    figure;contourf(wrs,tau,abs(direct),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Magnitude';'of Direct Fifth Order'});
    
    figure;contourf(wrs,tau,real(cascade_3),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Real Component';'Third Order Cascades'});
    
    figure;contourf(wrs,tau,imag(cascade_3),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Imaginary Component';'Third Order Cascades'});
    
    figure;contourf(wrs,tau,abs(cascade_3),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Magnitude';'of Third Order Cascades'});

    figure;contourf(wrs,tau,real(cascade_3./direct),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Real Component';'of Ratio E^{(3)}_{cascade}/E^{(5)}_{direct}'});
    
    figure;contourf(wrs,tau,imag(cascade_3./direct),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Imaginary Component';'of Ratio E^{(3)}_{cascade}/E^{(5)}_{direct}'});
    
    figure;contourf(wrs,tau,abs(cascade_3./direct),50,'edgecolor','none');
    axis square;xlabel('\omega (cm^{-1})');ylabel('\tau_1');xlim([0 max(wrs)]);
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title({'Magnitude of Ratio';'E^{(3)}_{cascade}/E^{(5)}_{direct}'});

end