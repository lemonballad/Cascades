% Load data from mat file
clear all
load FSRS
FSRS_Res_Cas=rcas;
FSRS_Res_Dir=rdir;
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.84 0.9]);%[0 4/75 1 0.85]);
maxr=max(max(max(abs(FSRS_Res_Ratio))));
y=-0.04;x=0.04;
% Axes positions
gcp_t_y=0.615;gcp_b_y=0.15;
gcp_l_x=0.05;gcp_c_x=0.32;gcp_r_x=0.59;
gcw=0.3;gch=0.3;
gcpos=[gcp_l_x gcp_t_y gcw gch;gcp_l_x gcp_b_y gcw gch];
gdpos=[gcp_c_x gcp_t_y gcw gch;gcp_c_x gcp_b_y gcw gch];
grpos=[gcp_r_x gcp_t_y gcw gch;gcp_r_x gcp_b_y gcw gch];
% Color bar positions
cp_t_y=0.925;cp_b_y=0.46;
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
for iv=1:2
    maxs=max(max(max(abs(FSRS_Res_Cas(:,:,iv)))),...
        max(max(abs(FSRS_Res_Dir(:,:,iv)))));
    gc(iv)=subplot(2,3,1+3*(iv-1));contourf(FSRS_Res_dws,FSRS_Res_disps,...
        abs(FSRS_Res_Cas(:,:,iv))/maxs,50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement, \it d';'Scaling Factor'});%caxis([0 maxr])
    colormap jet;
    cc(iv)=colorbar('Location','northoutside','fontsize',10,...
        'Position',ccpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gcpos(iv,:));
    axis square;
    tx=annotation(gcf,'textbox',tpos(1+3*(iv-1),:),...
        'String',tstr(1+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);

    gd(iv)=subplot(2,3,2+3*(iv-1));contourf(FSRS_Res_dws,FSRS_Res_disps,...
        abs(FSRS_Res_Dir(:,:,iv))/maxs,50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement, \it d';'Scaling Factor'});%caxis([0 maxr])
    colormap jet;cd(iv)=colorbar('Location','northoutside','fontsize',10,'Position',cdpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gdpos(iv,:));axis square;
    tx=annotation(gcf,'textbox',tpos(2+3*(iv-1),:),...
        'String',tstr(2+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);

    gr(iv)=subplot(2,3,3+3*(iv-1));contourf(FSRS_Res_dws,FSRS_Res_disps,...
        log10(FSRS_Res_Ratio(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement, \it d';'Scaling Factor'});%caxis([0 maxr])
    colormap jet;cr(iv)=colorbar('Location','northoutside','fontsize',10,'Position',crpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',grpos(iv,:));axis square;
    tx=annotation(gcf,'textbox',tpos(3+3*(iv-1),:),...
        'String',tstr(3+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);
end