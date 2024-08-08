% Load data from mat file
clear all
load FSRS3
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 0.8 0.9]);%[0 4/75 1 0.85]);
dw=1/FSRS_Res_nt;c=3e-5;
maxr=max(max(max(abs(FSRS_Res_Ratio))));
y=-0.04;x=0.04;
% Axes positions
gcp_t_y=0.585;gcp_b_y=0.10;
gcp_l_x=0.05;gcp_c_x=0.32;gcp_r_x=0.59;
gcw=0.3;gch=0.3;
gcpos=[gcp_l_x gcp_t_y gcw gch;gcp_l_x gcp_b_y gcw gch];
gdpos=[gcp_c_x gcp_t_y gcw gch;gcp_c_x gcp_b_y gcw gch];
grpos=[gcp_r_x gcp_t_y gcw gch;gcp_r_x gcp_b_y gcw gch];
% Color bar positions
cp_t_y=0.895;cp_b_y=0.41;
cp_l_x=0.1;cp_c_x=0.37;cp_r_x=0.64;
cw=0.2;ch=0.03;
ccpos=[cp_l_x cp_t_y cw ch;cp_l_x cp_b_y cw ch];
cdpos=[cp_c_x cp_t_y cw ch;cp_c_x cp_b_y cw ch];
crpos=[cp_r_x cp_t_y cw ch;cp_r_x cp_b_y cw ch];
% Text
tstr={'a)';'b)';'c)';'d)';'e)';'f)'};
% Text positions
tp_t_y=0.82;tp_b_y=0.33;
tp_l_x=0.105;tp_c_x=0.375;tp_r_x=0.645;
tw=0.06;th=0.06;
tpos=[tp_l_x tp_t_y tw th; tp_c_x tp_t_y tw th; tp_r_x tp_t_y tw th;...
          tp_l_x tp_b_y tw th; tp_c_x tp_b_y tw th; tp_r_x tp_b_y tw th];
      wv=[861 1600];
for iv=1:2
    maxs=max(max(max(abs(FSRS_Res_Cascade(:,:,iv)))),...
        max(max(abs(FSRS_Res_Direct(:,:,iv)))));
%     dt=dts(iv);
    w=(-1/2:dw:(1/2-dw))/FSRS_Res_dts(iv)/c;w=w-w(1);
    t=0:FSRS_Res_dts(iv):(FSRS_Res_nt-1)*FSRS_Res_dts(iv);
    if iv==1;wrange=[-1500 1500];else wrange=[-2400 2400];end
    gc(iv)=subplot(2,3,1+3*(iv-1));contourf(FSRS_Res_dws,...
        FSRS_Res_disps,abs(rcas(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode';'Displacement \it d'});%caxis([0 maxr])
    colormap jet;cc(iv)=colorbar('Location','northoutside','fontsize',10,'Position',ccpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gcpos(iv,:));
    axis square;
    legend('CASCADES','fontsize',10,'fontweight','bold')
    tx=annotation(gcf,'textbox',tpos(1+3*(iv-1),:),...
        'String',tstr(1+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);
    
    gd(iv)=subplot(2,3,2+3*(iv-1));contourf(FSRS_Res_dws,...
        FSRS_Res_disps,abs(rdir(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    legend('DIRECT','fontsize',10,'fontweight','bold')
    ylabel({'Mode';'Displacement \it d'});%caxis([0 maxr])
    colormap jet;cd(iv)=colorbar('Location','northoutside','fontsize',10,'Position',cdpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',gdpos(iv,:));axis square;%title('Direct')
    annotation(gcf,'textbox',tpos(2+3*(iv-1),:),...
        'String',tstr(2+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    gr(iv)=subplot(2,3,3+3*(iv-1));contourf(FSRS_Res_dws,...
        FSRS_Res_disps,(FSRS_Res_Ratio(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    legend('RATIO','fontsize',10,'fontweight','bold')
    ylabel({'Mode';'Displacement \it d'});%caxis([0.2 1])
    colormap jet;cr(iv)=colorbar('Location','northoutside','fontsize',10,'Position',crpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold','Position',grpos(iv,:));axis square;%title('E^{3}_{CAS}:E^{5}_{DIRECT}')
    annotation(gcf,'textbox',tpos(3+3*(iv-1),:),...
        'String',tstr(3+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);
end
% w_ap=weg+dw;w_rp=w_ap;
% w_t=w_rp-wv;
 
% L_ap_p=1./(dw+1i*gamma_eg);
% L_ap_m=1./(-dw+1i*gamma_eg);
% L_rp_m=1./(-dw+1i*gamma_eg);
% L_t=1./(dw-wv+1i*gamma_eg);

% D=1/LAMBDA_ap;
% J=1./(-wv+1i*(gamma_vib-LAMBDA_rp)).*exp(1i*(-wv)*r2w*tau2-LAMBDA_rp*tau2*r2w);


