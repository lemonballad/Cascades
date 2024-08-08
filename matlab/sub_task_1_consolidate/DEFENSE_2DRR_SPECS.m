% Load data from mat file
clear all
load Full2dRR

Cascade=Full_2dRR_OffRes_Cascade+Full_2dRR_Res_Cascade;
Direct=Full_2dRR_Res_Direct;
Full_Signal=Cascade+Direct;

dw=1/Full_2dRR_nt;c=3e-5;

      
    maxs=max(max(max(abs(real(Full_Signal(:,:,1))))));
    w=(-1/2:dw:(1/2-dw))/Full_2dRR_dts/c;
    figure;contour(w,w,abs(real(Full_2dRR_OffRes_Cascade(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');
    axis square;
    figure;contour(w,w,abs(real(Full_2dRR_Res_Cascade(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');
    axis square;
    figure;contour(w,w,abs(real(Cascade(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');
    axis square;
    
    figure;contour(w,w,abs(real(Direct(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;

    figure;contour(w,w,abs(real(Full_Signal(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;

