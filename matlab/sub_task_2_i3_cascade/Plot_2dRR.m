maxr=max(max(max(abs(ratio))));
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 4/75 1 0.855238095238095]);
for iv=1:2
    maxs=max(max(max(abs(cascade2d(:,:,iv)))),...
        max(max(abs(direct2d(:,:,iv)))));
    gc(iv)=subplot(2,3,1+3*(iv-1));contour(w,w,abs(cascade2d(:,:,iv))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');xlim([-2000 2000]);
    ylabel('\omega_2/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;colorbar('Location','northoutside','fontsize',16);
    set(gca,'linewidth',2,'fontsize',16);
    axis square;title('CASCADES')
    subplot(2,3,2+3*(iv-1));contour(w,w,abs(direct2d(:,:,iv))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');xlim([-2000 2000]);
    ylabel('\omega_2/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;colorbar('Location','northoutside','fontsize',16);
    set(gca,'linewidth',2,'fontsize',16);axis square;title('Direct')
    subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',16);
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end