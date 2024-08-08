% Load data from mat file
clear all
load 2dRR2d
load FSRS2d
set(0,'units','normalized');
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.3 0.9]);

% Axes positions
p_2dRR_y=0.49;p_FSRS_y=0.0;
p_2dRR_x=0.2;p_FSRS_x=0.2;
w=0.55;h=0.55;
pos=[p_2dRR_x p_2dRR_y w h;p_FSRS_x p_FSRS_y w h];

% Color bar positions
cp_2dRR_y=0.93;cp_FSRS_y=0.44;
cp_2dRR_x=0.2;cp_FSRS_x=0.2;
cw=0.5506;ch=0.03;
cpos=[cp_2dRR_x cp_2dRR_y cw ch;cp_FSRS_x cp_FSRS_y cw ch];

% Text
tstr={'a)';'b)'};
% Text positions
tp_2dRR_y=0.85;tp_FSRS_y=0.36;
tp_2dRR_x=0.2;tp_FSRS_x=0.2;
tw=0.06;th=0.06;
tpos=[tp_2dRR_x tp_2dRR_y tw th; tp_FSRS_x tp_FSRS_y tw th];

% Plot
subplot(2,1,1);contourf(FSRS_Res_2d_dws,FSRS_Res_2d_dvibs,(FSRS_Res_2d_Ratio(:,:)),50,'edgecolor','none');
xlabel('(\omega_L - \omega_{eg})/2\pic (cm^{-1})');
ylabel({'\omega_{vib}/2\pic (cm^{-1})'});%caxis([0 maxr])
colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cpos(1,:));
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',pos(1,:));axis square;
t1=annotation(gcf,'textbox',tpos(1,:),...
    'String',tstr(1),'EdgeColor','none',...
    'fontweight','bold','fontsize',10,'Color',[1 1 1]);

subplot(2,1,2);contourf(Res_2dRR_2d_dws,Res_2dRR_2d_dvibs,Res_2dRR_2d_Ratio(:,:),50,'edgecolor','none');
xlabel('(\omega_L - \omega_{eg})/2\pic (cm^{-1})');
ylabel({'\omega_{vib}/2\pic (cm^{-1})'});%caxis([0 maxr])
colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cpos(2,:));
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',pos(2,:));axis square;
t2=annotation(gcf,'textbox',tpos(2,:),...
    'String',tstr(2),'EdgeColor','none',...
    'fontweight','bold','fontsize',10,'Color',[1 1 1]);

