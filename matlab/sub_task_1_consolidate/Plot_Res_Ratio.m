% Load data from mat file
clear all
load 2dRR2d
load FSRS2d
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.84 0.9]);%[0 4/75 1 0.85]);

% Axes positions
gcp_t_y=0.61;gcp_b_y=0.1;
gcp_l_x=0.05;gcp_c_x=0.32;gcp_r_x=0.59;
gcw=0.3;gch=0.3;
gcpos=[gcp_l_x gcp_t_y gcw gch;gcp_l_x gcp_b_y gcw gch];
gdpos=[gcp_c_x gcp_t_y gcw gch;gcp_c_x gcp_b_y gcw gch];
grpos=[gcp_r_x gcp_t_y gcw gch;gcp_r_x gcp_b_y gcw gch];
% Color bar positions
cp_t_y=0.92;cp_b_y=0.41;
cp_l_x=0.1;cp_c_x=0.37;cp_r_x=0.64;
cw=0.2;ch=0.03;
ccpos=[cp_l_x cp_t_y cw ch;cp_l_x cp_b_y cw ch];
cdpos=[cp_c_x cp_t_y cw ch;cp_c_x cp_b_y cw ch];
crpos=[cp_r_x cp_t_y cw ch;cp_r_x cp_b_y cw ch];
% Text
tstr={'a)';'b)';'c)';'d)';'e)';'f)'};
% Text positions
tp_t_y=0.85;tp_b_y=0.38;
tp_l_x=0.105;tp_c_x=0.375;tp_r_x=0.645;
tw=0.06;th=0.06;
tpos=[tp_l_x tp_t_y tw th; tp_c_x tp_t_y tw th; tp_r_x tp_t_y tw th;...
          tp_l_x tp_b_y tw th; tp_c_x tp_b_y tw th; tp_r_x tp_b_y tw th];
%     maxs=max(max(max(abs(Res_2dRR_Cas(:,:,iv)))),...
%         max(max(abs(Res_2dRR_Dir(:,:,iv)))));
    subplot(2,3,1);contourf(FSRS_Res_2d_dws,FSRS_Res_2d_dvibs,...
        abs(FSRS_Res_2d_Dir),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'(\omega_t-\omega_RP)/2\pic (cm^{-1})'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',ccpos(1,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gcpos(1,:));
    axis square;
    annotation(gcf,'textbox',tpos(1,:),...
        'String',tstr(1),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);
    
    subplot(2,3,2);contourf(FSRS_Res_2d_dws,FSRS_Res_2d_dvibs,...
        abs(FSRS_Res_2d_Cas),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'(\omega_t-\omega_RP)/2\pic (cm^{-1})'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cdpos(1,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gdpos(1,:));axis square;
    annotation(gcf,'textbox',tpos(2,:),...
        'String',tstr(2),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    subplot(2,3,3);contourf(FSRS_Res_2d_dws,FSRS_Res_2d_dvibs,(FSRS_Res_2d_Ratio),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'(\omega_t-\omega_RP)/2\pic (cm^{-1})'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',crpos(1,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',grpos(1,:));axis square;
    annotation(gcf,'textbox',tpos(3,:),...
        'String',tstr(3),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);
    
    subplot(2,3,4);contourf(Res_2dRR_2d_dws,Res_2dRR_2d_dvibs,...
        abs(Res_2dRR_2d_Dir),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'(\omega_t-\omega_RP)/2\pic (cm^{-1})'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',ccpos(2,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gcpos(2,:));
    axis square;
    annotation(gcf,'textbox',tpos(2,:),...
        'String',tstr(4),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);
    
    subplot(2,3,5);contourf(Res_2dRR_2d_dws,Res_2dRR_2d_dvibs,...
        abs(Res_2dRR_2d_Cas),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'(\omega_t-\omega_RP)/2\pic (cm^{-1})'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cdpos(2,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gdpos(1,:));axis square;
    annotation(gcf,'textbox',tpos(2,:),...
        'String',tstr(5),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    subplot(2,3,6);contourf(Res_2dRR_2d_dws,Res_2dRR_2d_dvibs,(Res_2dRR_2d_Ratio),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'(\omega_t-\omega_RP)/2\pic (cm^{-1})'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',crpos(2,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',grpos(1,:));axis square;
    annotation(gcf,'textbox',tpos(2,:),...
        'String',tstr(6),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);

    
