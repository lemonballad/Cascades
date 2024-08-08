% Load data from mat file
clear all
load 2dOffRR2d
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.3 0.9]);

% Axes positions
p_2dRR_t_y=0.49;p_2dRR_b_y=0.0;
p_2dRR_t_x=0.21;p_2dRR_b_x=0.21;
w=0.55;h=0.55;
pos=[p_2dRR_t_x p_2dRR_t_y w h;p_2dRR_b_x p_2dRR_b_y w h];

% Color bar positions
cp_2dRR_y=0.935;cp_2dRR_b_y=0.445;
cp_2dRR_x=0.21;cp_2dRR_b_x=0.21;
cw=0.5506;ch=0.03;
cpos=[cp_2dRR_x cp_2dRR_y cw ch;cp_2dRR_b_x cp_2dRR_b_y cw ch];

% Text
tstr={'a)';'b)'};
% Text positions
tp_2dRR_y=0.85;tp_2dRR_b_y=0.36;
tp_2dRR_x=0.21;tp_2dRR_b_x=0.21;
tw=0.06;th=0.06;
tpos=[tp_2dRR_x tp_2dRR_y tw th; tp_2dRR_b_x tp_2dRR_b_y tw th];

% Plot
subplot(2,1,1);contourf(OffRes_2dRR_2d_dws,OffRes_2dRR_2d_dvibs,OffRes_2dRR_2d_Ratio(:,:,1),50,'edgecolor','none');
xlabel('(\omega_L - \omega_{eg})/2\pic (cm^{-1})');
ylabel({'\omega_{vib}/2\pic (cm^{-1})'});%caxis([0 maxr])
colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cpos(1,:));
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',pos(1,:));axis square;
annotation(gcf,'textbox',tpos(1,:),...
    'String',tstr(1),'EdgeColor','none',...
    'fontweight','bold','fontsize',10,'Color',[1 1 1]);

subplot(2,1,2);contourf(OffRes_2dRR_2d_dws,OffRes_2dRR_2d_dvibs,OffRes_2dRR_2d_Ratio(:,:,2),50,'edgecolor','none');
xlabel('(\omega_L - \omega_{eg})/2\pic (cm^{-1})');
ylabel({'\omega_{vib}/2\pic (cm^{-1})'});%caxis([0 maxr])
colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cpos(2,:));
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',pos(2,:));axis square;
annotation(gcf,'textbox',tpos(2,:),...
    'String',tstr(2),'EdgeColor','none',...
    'fontweight','bold','fontsize',10,'Color',[1 1 1]);