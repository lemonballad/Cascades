figure;subplot(2,2,1);contourf(w_aps,disps,abs(cas(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(w_aps,disps,abs(cas2(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(w_aps,disps,abs(cas(:,:,1))./abs(cas2(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(w_aps,disps,abs(cas(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(w_aps,disps,abs(cas2(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(w_aps,disps,abs(cas(:,:,2))./abs(cas2(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;


figure;subplot(2,2,1);contourf(w_aps,disps,abs(dir(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(w_aps,disps,abs(dir2(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(w_aps,disps,abs(dir(:,:,1))./abs(dir2(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(w_aps,disps,abs(dir(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(w_aps,disps,abs(dir2(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(w_aps,disps,abs(dir(:,:,2))./abs(dir2(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;


figure;subplot(2,2,1);contourf(w_aps,disps,abs(ratio(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(w_aps,disps,abs(ratio2(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(w_aps,disps,abs(ratio(:,:,1))./abs(ratio2(:,:,1)),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(w_aps,disps,abs(ratio(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(w_aps,disps,abs(ratio2(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(w_aps,disps,abs(ratio(:,:,2))./abs(ratio2(:,:,2)),50,'edgecolor','none');colormap('jet');colorbar;

