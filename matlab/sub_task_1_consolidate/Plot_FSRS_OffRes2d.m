% Load data from mat file
clear all
load FSRSOff2d
set(0,'units','normalized');
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.3 0.9]);

% Axes positions
p_FSRS_t_y=0.49;p_FSRS_b_y=0.0;
p_FSRS_t_x=0.2;p_FSRS_b_x=0.2;
w=0.55;h=0.55;
pos=[p_FSRS_t_x p_FSRS_t_y w h;p_FSRS_b_x p_FSRS_b_y w h];

% Color bar positions
cp_FSRS_t_y=0.935;cp_FSRS_y=0.445;
cp_FSRS_t_x=0.2;cp_FSRS_x=0.2;
cw=0.5506;ch=0.03;
cpos=[cp_FSRS_t_x cp_FSRS_t_y cw ch;cp_FSRS_x cp_FSRS_y cw ch];

% Text
tstr={'a)';'b)'};
% Text positions
tp_FSRS_t_y=0.85;tp_FSRS_y=0.36;
tp_FSRS_t_x=0.2;tp_FSRS_x=0.2;
tw=0.06;th=0.06;
tpos=[tp_FSRS_t_x tp_FSRS_t_y tw th; tp_FSRS_x tp_FSRS_y tw th];

% Plot
subplot(2,1,1);contourf(FSRS_OffRes_2d_dws,...
    FSRS_OffRes_2d_dvibs,log10(abs(FSRS_OffRes_2d_Ratio(:,:,1))),50,'edgecolor','none');
xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
ylabel({'\omega_{vib}/2\pic (cm^{-1})'});
colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cpos(1,:));
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',pos(1,:));axis square;
annotation(gcf,'textbox',tpos(1,:),...
    'String',tstr(1),'EdgeColor','none',...
    'fontweight','bold','fontsize',10,'Color',[1 1 1]);

subplot(2,1,2);contourf(FSRS_OffRes_2d_dws,...
    FSRS_OffRes_2d_dvibs,log10(abs(FSRS_OffRes_2d_Ratio(:,:,2))),50,'edgecolor','none');
xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
ylabel({'\omega_{vib}/2\pic (cm^{-1})'});
colormap jet;colorbar('Location','northoutside','fontsize',10,'Position',cpos(2,:));
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',pos(2,:));axis square;
annotation(gcf,'textbox',tpos(2,:),...
    'String',tstr(1),'EdgeColor','none',...
    'fontweight','bold','fontsize',10,'Color',[1 1 1]);