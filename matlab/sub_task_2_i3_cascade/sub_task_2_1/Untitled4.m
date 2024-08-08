figure
for iv=1:2    
subplot(1,2,iv);contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');
    colormap jet;colorbar('Location','northoutside','fontsize',16);%caxis([0 0.2])
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end


maxr=max(max(max(abs(ratio))));
    figure;
for iv=1:2
    maxc=max(max(abs(cascade2d(:,:,iv))));
    maxd=max(max(abs(direct2d(:,:,iv))));
    maxcd=max(maxc,maxd);
    subplot(2,3,1+3*(iv-1));contour(w,w,abs(cascade2d(:,:,iv))/maxcd,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',12);caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('CASCADES')
    subplot(2,3,2+3*(iv-1));contour(w,w,abs(direct2d(:,:,iv))/maxcd,50);
    xlabel('\omega_1/2\pic (cm^{-1})');
    %ylabel('\omega_2/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',12);caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12,'ytick',[]);axis square;%title('Direct')
    subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode';'Displacement'});
    colormap jet;colorbar('Location','northoutside','fontsize',12);caxis([0 maxr]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('E^{3}_{CAS}:E^{5}_{DIRECT}');
end




rTC=prefactor*abs(cascade_3)./abs(direct)/3e10;
figure
for iv=1:2    
subplot(1,2,iv);contourf(w_aps,disps,rTC(:,:,iv,1),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');
    colormap jet;colorbar('Location','northoutside','fontsize',16);%caxis([0 0.2])
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end



rTC=prefactor*abs(cascade_3)./abs(direct)/3e10;
figure
for iv=1:2
subplot(1,2,iv);contour(w_aps,disps,rTC(:,:,iv,85),50);
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');xlim([-1000 1000])
    colormap jet;colorbar('Location','northoutside','fontsize',16);%caxis([0 0.2])
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end




    figure;contour(wvibs,disps,ratio_AM,35);%     axis square;
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');
    colormap jet;colorbar('Location','northoutside','fontsize',16)
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')

