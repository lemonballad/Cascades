figure;subplot(2,2,1);contourf(dws,disps,rcasG(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rcasG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,rcasG(:,:,1)./rcasG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(dws,disps,rcasG(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rcasG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,rcasG(:,:,2)./rcasG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;

figure;subplot(2,2,1);contourf(dws,disps,ratioG(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,ratioG(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;

figure;subplot(2,2,1);contourf(dws,disps,ratioG(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,ratioG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,ratioG(:,:,1)./ratioG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(dws,disps,ratioG(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,ratioG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,ratioG(:,:,2)./ratioG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;


figure;subplot(2,2,1);contourf(dws,disps,rdirG(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rdirG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,rdirG(:,:,1)./rdirG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(dws,disps,rdirG(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rdirG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,rdirG(:,:,2)./rdirG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;



figure;subplot(2,2,1);contourf(dws,disps,rcasG(:,:,1)./rdirG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rcasG2(:,:,1)./rdirG(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,ratioG(:,:,1)./ratioG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(dws,disps,rcasG(:,:,2)./rdirG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rcasG2(:,:,2)./rdirG(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,ratioG(:,:,2)./ratioG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;


figure;plot(rcasG(4,:,1),rcasG2(4,:,1),'k-','linewidth',2);
figure;plot(rcasG(4,:,2),rcasG2(4,:,2),'k-','linewidth',2);
figure;plot(rdirG(4,:,1),rdirG2(4,:,1),'k-','linewidth',2);
figure;plot(rdirG(4,:,2),rdirG2(4,:,2),'k-','linewidth',2);
figure;plot(ratioG(4,:,1),ratioG2(4,:,1),'k-','linewidth',2);
figure;plot(ratioG(4,:,2),ratioG2(4,:,2),'k-','linewidth',2);


figure;plot(permute(F(1,4,:,1),[3 1 2 4]),permute(F2(1,4,:,1),[3 1 2 4]),'k-','linewidth',2);
figure;plot(permute(F(2,4,:,1),[3 1 2 4]),permute(F2(2,4,:,1),[3 1 2 4]),'k-','linewidth',2);
figure;plot(permute(F(3,4,:,1),[3 1 2 4]),permute(F2(3,4,:,1),[3 1 2 4]),'k-','linewidth',2);
figure;plot(permute(F(4,4,:,1),[3 1 2 4]),permute(F2(4,4,:,1),[3 1 2 4]),'k-','linewidth',2);
figure;plot(permute(F(5,4,:,1),[3 1 2 4]),permute(F2(5,4,:,1),[3 1 2 4]),'k-','linewidth',2);
