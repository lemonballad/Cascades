% Load data from mat file
clear all
load FSRS
load 2dRR

id=16;find(FSRS_Res_disps==1);
figure;plot(FSRS_Res_dws,(FSRS_Res_Ratio(id,:,1)),'r-','Linewidth',2);hold on
plot(FSRS_Res_dws,(FSRS_Res_Ratio(id,:,2)),'b-','Linewidth',2);
xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

iw=find(FSRS_Res_dws==0);
figure;plot(FSRS_Res_disps,(FSRS_Res_Ratio(:,iw,1)),'r-','Linewidth',2);hold on
plot(FSRS_Res_disps,(FSRS_Res_Ratio(:,iw,2)),'b-','Linewidth',2);xlim([0.1 inf])
xlabel({'Mode Displacement, \it d';'Scaling Factor'});ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

id=16;find(Res_2dRR_disps==1);
figure;plot(Res_2dRR_dws,(Res_2dRR_Ratio(id,:,1)),'r-','Linewidth',2);hold on
plot(Res_2dRR_dws,(Res_2dRR_Ratio(id,:,2)),'b-','Linewidth',2);
xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

iw=find(Res_2dRR_dws==0);
figure;plot(Res_2dRR_disps,(Res_2dRR_Ratio(:,iw,1)),'r-','Linewidth',2);hold on
plot(Res_2dRR_disps,(Res_2dRR_Ratio(:,iw,2)),'b-','Linewidth',2);xlim([0.1 inf])
xlabel({'Mode Displacement, \it d';'Scaling Factor'});ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

