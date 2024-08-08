nt=100000;
dt=20;
t=0:dt:(nt-1)*dt;
dw=1/nt;
w=(-1/2:dw:(1/2-dw))/dt/c;
A=exp(-1i*500*2*pi*c*t)*nt*2*pi*dt/(500*c);
A=fftshift(fft(A));A=abs(A);
B=1./(w-500);B=abs(B);
figure;plot(w,A,w,B);



figure;subplot(2,2,1);contourf(dws,disps,rcasG(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rcasG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,rcasG(:,:,1)./rcasG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(dws,disps,rcasG(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rcasG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,rcasG(:,:,2)./rcasG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;


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
subplot(2,2,3);contourf(dws,disps,ratioG(:,:,2)./rdirG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;



figure;subplot(2,2,1);contourf(dws,disps,rcasG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rdirG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,ratioG2(:,:,1),50,'edgecolor','none');colormap('jet');colorbar;
figure;subplot(2,2,1);contourf(dws,disps,rcasG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,2);contourf(dws,disps,rdirG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;
subplot(2,2,3);contourf(dws,disps,ratioG2(:,:,2),50,'edgecolor','none');colormap('jet');colorbar;