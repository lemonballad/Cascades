load FSRS2.m;
maxr=max(max(max(abs(ratio))));
% gcf=figure;
% set(gcf,'Units','Normalized','Position',[0 4/75 1 0.855238095238095]);
for iv=1:2
    maxs=max(max(max(abs(cascadeFSRS(:,:,iv)))),...
        max(max(abs(directFSRS(:,:,iv)))));
    subplot(2,3,1+3*(iv-1));contour(t/1e3,w,abs(cascadeFSRS(:,:,iv))'/maxs,50);
    xlabel('\tau_1 (ps)');xlim([0 7]);
    ylabel('\omega/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10);
    axis square;title('CASCADES')
    subplot(2,3,2+3*(iv-1));contour(t/1e3,w,abs(directFSRS(:,:,iv))'/maxs,50);
    xlabel('\tau_1 (ps)');xlim([0 7]);
    ylabel('\omega/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10);axis square;title('Direct')
    subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode';'Displacement'});caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',10);
    set(gca,'linewidth',2,'fontsize',10);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end