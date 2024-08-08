clear all
load 2dOffRR2d
load FSRSOff2d

iw_FSRS=find(FSRS_OffRes_2d_dws==0);
iw_2dRR=find(OffRes_2dRR_2d_dws==0);


figure;plot(FSRS_OffRes_2d_dvibs,FSRS_OffRes_2d_Ratio(:,iw_FSRS,1),'r-',...
    FSRS_OffRes_2d_dvibs,FSRS_OffRes_2d_Ratio(:,iw_FSRS,2),'b-',...
    'linewidth',2)
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
l=legend('861 cm^{-1}','1328 cm^{-1}');l.Position=[0.545 0.7992 0.215 0.0952];
xlabel({'\omega_{vib}/2\pic (cm^{-1})'});
ylabel({'|E^{(3)}_{Solute-Solvent,Cas}|/|E^{(5)}_{FSRS}|'});
xlim([500 FSRS_OffRes_2d_dvibs(end)]);

max_FSRS=max(max(max(FSRS_OffRes_2d_Ratio)));
p = get(gca, 'Position');
h_861 = axes('Parent', gcf, 'Position', [p(1)+.396 p(2)+.4 p(3)-.55 p(4)-.55]);
contourf(h_861,FSRS_OffRes_2d_dws,FSRS_OffRes_2d_dvibs,...
    FSRS_OffRes_2d_Ratio(:,:,1),50,'edgecolor','none');
y=ylabel({'\omega_{vib}/2\pic (cm^{-1})'},'position',[-4.5085 0.7 0.001]*1e3);
colormap jet;caxis([0 max_FSRS]);
c=colorbar('fontsize',8,'position',[0.74 0.254 0.02 0.523]);
set(h_861,'linewidth',2,'fontsize',8,'fontweight','bold',...
    'ytick',861,'xticklabel',{[]});
 axis square;
h_1328 = axes('Parent', gcf, 'Position', [p(1)+.396 p(2)+.14 p(3)-.55 p(4)-.55]);
contourf(FSRS_OffRes_2d_dws,FSRS_OffRes_2d_dvibs,...
    FSRS_OffRes_2d_Ratio(:,:,2),50,'edgecolor','none');
xlabel('(\omega_{AP} - \omega_{eg})/2\pic (cm^{-1})');
colormap jet;caxis([0 max_FSRS]);
set(h_1328,'linewidth',2,'fontsize',8,'fontweight','bold',...
    'ytick',1328);
 axis square;


% 2dRR
figure;plot(OffRes_2dRR_2d_dvibs,OffRes_2dRR_2d_Ratio(:,iw_2dRR,1),'r-',...
    OffRes_2dRR_2d_dvibs,OffRes_2dRR_2d_Ratio(:,iw_2dRR,2),'b-',...
    'linewidth',2)
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
l=legend('861 cm^{-1}','1328 cm^{-1}');l.Position=[0.545 0.7992 0.215 0.0952];
xlabel({'\omega_{vib}/2\pic (cm^{-1})'});
ylabel({'|E^{(3)}_{Solute-Solvent,Cas}|/|E^{(5)}_{2DRR}|'});
xlim([500 OffRes_2dRR_2d_dvibs(end)]);

max_2dRR=max(max(max(OffRes_2dRR_2d_Ratio)));
p = get(gca, 'Position');
h_861 = axes('Parent', gcf, 'Position', [p(1)+.396 p(2)+.4 p(3)-.55 p(4)-.55]);
contourf(h_861,OffRes_2dRR_2d_dws,OffRes_2dRR_2d_dvibs,...
    OffRes_2dRR_2d_Ratio(:,:,1),50,'edgecolor','none');
y=ylabel({'\omega_{vib}/2\pic (cm^{-1})'},'position',[-4.5085 0.7 0.001]*1e3);
colormap jet;caxis([0 max_2dRR]);
c=colorbar('fontsize',8,'position',[0.74 0.254 0.02 0.523]);
set(h_861,'linewidth',2,'fontsize',8,'fontweight','bold',...
    'ytick',861,'xticklabel',{[]});
 axis square;
h_1328 = axes('Parent', gcf, 'Position', [p(1)+.396 p(2)+.14 p(3)-.55 p(4)-.55]);
contourf(OffRes_2dRR_2d_dws,OffRes_2dRR_2d_dvibs,...
    OffRes_2dRR_2d_Ratio(:,:,2),50,'edgecolor','none');
xlabel('(\omega_L - \omega_{eg})/2\pic (cm^{-1})');
colormap jet;caxis([0 max_2dRR]);
set(h_1328,'linewidth',2,'fontsize',8,'fontweight','bold',...
    'ytick',1328);
 axis square;
% annotation(gcf,'textbox',tpos(1,:),...
%     'String',tstr(1),'EdgeColor','none',...
%     'fontweight','bold','fontsize',10,'Color',[1 1 1]);
