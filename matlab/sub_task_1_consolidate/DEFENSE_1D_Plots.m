% Load data from mat file
clear all
load FSRSOff
load 2dOffRR
load FSRS
load 2dRR

id=find(FSRS_OffRes_disps==1);
figure;plot(FSRS_OffRes_dws,(FSRS_OffRes_Ratio(id,:,1)+...
    FSRS_Res_Ratio(id,:,1)),'r-','Linewidth',2);hold on
plot(FSRS_OffRes_dws,(FSRS_OffRes_Ratio(id,:,2)+...
    FSRS_Res_Ratio(id,:,2)),'b-','Linewidth',2);
xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

iw=find(FSRS_OffRes_dws==0);
figure;plot(FSRS_OffRes_disps,(FSRS_OffRes_Ratio(:,iw,1)+...
    FSRS_Res_Ratio(:,iw,1)),'r-','Linewidth',2);hold on
plot(FSRS_OffRes_disps,(FSRS_OffRes_Ratio(:,iw,2)+...
    FSRS_Res_Ratio(:,iw,2)),'b-','Linewidth',2);xlim([0.1 inf])
xlabel({'Mode Displacement, \it d';'Scaling Factor'});ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

id=find(OffRes_2dRR_disps==1);
figure;plot(OffRes_2dRR_dws,(OffRes_2dRR_Ratio(id,:,1)+...
    Res_2dRR_Ratio(id,:,1)),'r-','Linewidth',2);hold on
plot(OffRes_2dRR_dws,(OffRes_2dRR_Ratio(id,:,2)+...
    Res_2dRR_Ratio(id,:,2)),'b-','Linewidth',2);
xlabel('(\omega_{AP}-\omega_{eg})/2\pic (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

iw=find(OffRes_2dRR_dws==0);
figure;plot(OffRes_2dRR_disps,(OffRes_2dRR_Ratio(:,iw,1)+...
    Res_2dRR_Ratio(:,iw,1)),'r-','Linewidth',2);hold on
plot(OffRes_2dRR_disps,(OffRes_2dRR_Ratio(:,iw,2)+...
    Res_2dRR_Ratio(:,iw,2)),'b-','Linewidth',2);xlim([0.1 inf])
xlabel({'Mode Displacement, \it d';'Scaling Factor'});ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
hold off

