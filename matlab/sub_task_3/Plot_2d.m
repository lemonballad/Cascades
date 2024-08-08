load 2dRR
maxr=max(max(max(abs(ratio))));
gcpos=[0.05 0.59 0.33 0.33; 0.05 0.09 0.33 0.33];
gdpos=[0.33 0.59 0.33 0.33; 0.33 0.09 0.33 0.33];
grpos=[0.61 0.59 0.33 0.33; 0.61 0.09 0.33 0.33];
ccpos=[0.126 0.925 0.178 0.03;0.126 0.425 0.178 0.03];
cdpos=[0.406 0.925 0.178 0.03;0.406 0.425 0.178 0.03];
crpos=[0.686 0.925 0.178 0.03;0.686 0.425 0.178 0.03];
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 4/75 1 0.85]);
for iv=1:2
    maxs=max(max(max(prefactor(iv)/3e10*abs(cascade2d(:,:,iv)))),...
        max(max(abs(direct2d(:,:,iv)))));
    gc(iv)=subplot(2,3,1+3*(iv-1));contour(w,w,prefactor(iv)/3e10*abs(cascade2d(:,:,iv))'/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');xlim([-2000 2000]);
    ylabel('\omega_2/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;cc(iv)=colorbar('Location','northoutside','fontsize',10,'Position',ccpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'Position',gcpos(iv,:));
    axis square;%title('CASCADES')
    gd(iv)=subplot(2,3,2+3*(iv-1));contour(w,w,abs(direct2d(:,:,iv))'/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');xlim([-2000 2000]);
    ylabel('\omega_2/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;cd(iv)=colorbar('Location','northoutside','fontsize',10,'Position',cdpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'Position',gdpos(iv,:));axis square;%title('Direct')
    gr(iv)=subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode';'Displacement'});caxis([0 maxr])
    colormap jet;cr(iv)=colorbar('Location','northoutside','fontsize',10,'Position',crpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'Position',grpos(iv,:));axis square;%title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end