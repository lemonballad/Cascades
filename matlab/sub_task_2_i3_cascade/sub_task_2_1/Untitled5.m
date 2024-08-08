% maxr=max(max(max(abs(ratio))));
% t=0:dt:nt*dt-dt;
% dw=1/nt;
% w=(-1/2:dw:(1/2-dw))/dt/c;
% [~,iomega]=min(abs(w-wvib));
for iv=1:2
figure;
    maxc=max(max(abs(cascade_3(:,:,iv))));
    maxd=max(max(abs(direct(:,:,iv))));
    maxcd=max(maxc,maxd);
    subplot(1,2,1);contour(t/1e3,w,prefactor(iv)/3e10*abs(cascade_3(:,:,iv))'/maxcd,50);
    xlabel('\tau_1 (ps)');xlim([0 7]);
    ylabel('\omega/2\pic (cm^{-1})');ylim([0 2000]);
    colormap jet;colorbar('Location','northoutside','fontsize',12);%caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('CASCADES')
    subplot(1,2,2);contour(t/1e3,w,abs(direct(:,:,iv))'/maxcd,50);
    xlabel('\tau_1 (ps)');xlim([0 7]);
    ylabel('\omega/2\pic (cm^{-1})');ylim([0 2000]);
    colormap jet;colorbar('Location','northoutside','fontsize',12);%caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('Direct')
%     subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
%     xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
%     ylabel({'Mode';'Displacement'});
%     colormap jet;colorbar('Location','northoutside','fontsize',12);caxis([0 maxr]);
%     set(gca,'linewidth',2,'fontsize',12);axis square;%title('E^{3}_{CAS}:E^{5}_{DIRECT}');
end


fc=fftshift(fft(cascade_3,[],1),1);
fd=fftshift(fft(direct,[],1),1);
for iv=1:2
figure;
    maxc=max(max(abs(fc(:,:,iv))));
    maxd=max(max(abs(fd(:,:,iv))));
    maxcd=max(maxc,maxd);
    subplot(1,2,1);contourf(w-w(end)/2,w,abs(fc(:,:,iv))'/maxcd,50,'edgecolor','none');
    xlabel('\tau_1 (ps)');
    ylabel('\omega/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',12);%caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('CASCADES')
    subplot(1,2,2);contourf(w-w(end)/2,w,abs(fd(:,:,iv))'/maxcd,50,'edgecolor','none');
    xlabel('\tau_1 (ps)');
    ylabel('\omega/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',12);%caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('Direct')
%     subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
%     xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
%     ylabel({'Mode';'Displacement'});
%     colormap jet;colorbar('Location','northoutside','fontsize',12);caxis([0 maxr]);
%     set(gca,'linewidth',2,'fontsize',12);axis square;%title('E^{3}_{CAS}:E^{5}_{DIRECT}');
end




figure; plot(w,abs(direct(1,:,2)));
figure; plot(t,abs(direct(:,427,1)));



figure;
    maxc=max(max(abs(seq))); maxd=max(max(abs(direct))); maxcd=max(maxc,maxd);
    subplot(1,2,1);contour(tau/1e3,wrs,abs(seq)'/maxcd,50);
    xlabel('\tau_1 (ps)');    ylabel('\omega/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',12);%caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('CASCADES')
    subplot(1,2,2);contour(tau/1e3,wrs,abs(direct)'/maxcd,50);
    xlabel('\tau_1 (ps)');    ylabel('\omega/2\pic (cm^{-1})');
    colormap jet;colorbar('Location','northoutside','fontsize',12);%caxis([0 1]);
    set(gca,'linewidth',2,'fontsize',12);axis square;%title('Direct')


    rTC=prefactor.*abs(cascade_3(:,:,:,1))./abs(direct(:,:,:,1))/3e10;
    rTC(isnan(rTC))=0;
    figure;
    for iv=1:2
    subplot(1,2,iv);contourf(w_aps,disps,rTC(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel('Mode Displacement');%caxis([0 maxr])
    colormap jet;colorbar('Location','northoutside','fontsize',16);
    set(gca,'linewidth',2,'fontsize',16);axis square;title('E^{3}_{CAS}:E^{5}_{DIRECT}')
    end
    