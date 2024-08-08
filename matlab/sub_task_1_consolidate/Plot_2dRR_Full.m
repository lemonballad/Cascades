% Load data from mat file
clear all
load Full2dRR

Cascade=Full_2dRR_OffRes_Cascade+Full_2dRR_Res_Cascade;
Direct=Full_2dRR_Res_Direct;
Full_Signal=Cascade+Direct;

gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.85 0.9]);%[0 4/75 1 0.85]);
dw=1/Full_2dRR_nt;c=3e-5;

% Axes positions
gcp_y=0.3;
gcp_l_x=0.05;gcp_c_x=0.32;gcp_r_x=0.59;
gcw=0.3;gch=0.3;
gcpos=[gcp_l_x gcp_y gcw gch];
gdpos=[gcp_c_x gcp_y gcw gch];
grpos=[gcp_r_x gcp_y gcw gch];
% Color bar positions
cp_y=0.61;
cp_l_x=0.104;cp_c_x=0.374;cp_r_x=0.644;
cw=0.192;ch=0.03;
ccpos=[cp_l_x cp_y cw ch];
cdpos=[cp_c_x cp_y cw ch];
crpos=[cp_r_x cp_y cw ch];
% Text
tstr={'a)';'b)';'c)'};
% Text positions
tp_y=0.535;
tp_l_x=0.105;tp_c_x=0.375;tp_r_x=0.645;
tw=0.06;th=0.06;
tpos=[tp_l_x tp_y tw th; tp_c_x tp_y tw th; tp_r_x tp_y tw th];
      
    maxs=max(max(max(abs(real(Full_Signal(:,:,1))))));
%     dt=dts(iv);
    w=(-1/2:dw:(1/2-dw))/Full_2dRR_dts/c;
    subplot(1,3,1);contour(w,w,abs(real(Cascade(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',ccpos);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gcpos);
    axis square;
    annotation(gcf,'textbox',tpos(1,:),...
        'String',tstr(1),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);
    
    subplot(1,3,2);contour(w,w,abs(real(Direct(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cdpos);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gdpos);axis square;
    annotation(gcf,'textbox',tpos(2,:),...
        'String',tstr(2),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    subplot(1,3,3);contour(w,w,abs(real(Full_Signal(:,:,1)))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',crpos);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',grpos);axis square;
    annotation(gcf,'textbox',tpos(3,:),...
        'String',tstr(3),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);
