% Load data from mat file
clear all
load 2dRR
for iv=1:2
    figure
    contourf(Res_2dRR_dws,Res_2dRR_disps,...
        log10(Res_2dRR_Ratio(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement \it d';'Scaling Factor'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
end

% Load data from mat file
clear all
load FSRS
for iv=1:2
    figure
    contourf(FSRS_Res_dws,FSRS_Res_disps,...
        log10(FSRS_Res_Ratio(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement, \it d';'Scaling Factor'});%caxis([0 maxr])
    colormap jet;cr(iv)=colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
end

% Off Res

% Load data from mat file
clear all
load 2dOffRR
for iv=1:2

    figure
    contourf(OffRes_2dRR_dws,OffRes_2dRR_disps,...
        log10(OffRes_2dRR_Ratio(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement \it d';'Scaling Factor'});%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
end

% Load data from mat file
clear all
load FSRSOff
for iv=1:2
    figure
    contourf(FSRS_OffRes_dws,FSRS_OffRes_disps,...
        log10(FSRS_OffRes_Ratio(:,:,iv)),50,'edgecolor','none');
    xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode Displacement \it d';'Scaling Factor'});
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
end