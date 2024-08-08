figure
for iv=1:2    
subplot(1,2,iv);contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');
    colormap jet;colorbar('Location','northoutside','fontsize',16);%caxis([0 0.2])
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
end
