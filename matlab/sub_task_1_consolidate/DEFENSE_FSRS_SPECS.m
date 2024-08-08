% Load data from mat file
clear all
load FSRSFull
Cascade=1i*FSRS_Full_OffRes_Cascade+1*FSRS_Full_Res_Cascade;
Direct=-FSRS_Full_OffRes_Direct;
Full_Spec=Cascade+Direct;
dw=1/FSRS_Full_nt;c=3e-5;

    maxs=max(max(max(abs(real(Full_Spec(:,:,1))))),...
        max(max(abs(real(Direct(:,:,1))))));
    w=(-1/2:dw:(1/2-dw))/FSRS_Full_dts/c;w=w-w(1);
    t=0:FSRS_Full_dts:(FSRS_Full_nt-1)*FSRS_Full_dts;
%     figure;contour(t,w,abs(real(FSRS_Full_Res_Cascade(:,:,1)))/maxs,50);
%     xlabel('\tau (fs)');
%     ylabel('\omega/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',10);
%     set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');
%     axis square;
%     figure;contour(t,w,abs(real(1i*FSRS_Full_OffRes_Cascade(:,:,1)))/maxs,50);
%     xlabel('\tau (fs)');
%     ylabel('\omega/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',10);
%     set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');
%     axis square;
%     figure;contour(t,w,abs(real(Cascade(:,:,1)))/maxs,50);
%     xlabel('\tau (fs)');
%     ylabel('\omega/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',10);
%     set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');
%     axis square;
    
%     figure;contour(t,w,abs(real(Direct(:,:,1)))/maxs,50);
%     xlabel('\tau (fs)');
%     ylabel('\omega/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',10);
%     set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
% 
%     figure;contour(t,w,abs(real(Full_Spec(:,:,1)))/maxs,50);
%     xlabel('\tau (fs)');
%     ylabel('\omega/2\pic (cm^{-1})');
%     colormap jet;colorbar('Location','northoutside','fontsize',10);
%     set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;

    figure;plot(w,abs(real(Direct(:,end,1)))/maxs,'g-','Linewidth',2);
    xlabel('\omega/2\pic (cm^{-1})');xlim([0 3000]);ylabel('INTENSITY');
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;

    figure;plot(w,abs(real(Full_Spec(:,end,1)))/maxs,'g-','Linewidth',2);
    xlabel('\omega/2\pic (cm^{-1})');xlim([0 3000]);ylabel('INTENSITY');
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
