clear all
clc
load CH2Cl2_2DRR_Vib 
load FSRSCH2Cl2Vib

CH2Cl2_2DRR_Vib_wsolvs=0:50:2000;nw=length(CH2Cl2_2DRR_Vib_wsolvs);

figure;plot(CH2Cl2_2DRR_Vib_wsolvs,abs(rcas(:,1)+...
    rcas2(:,1))./abs(rdir(:,1)),'r-',CH2Cl2_2DRR_Vib_wsolvs,abs(CH2Cl2_2DRR_Vib_Cas(:,1)+...
    CH2Cl2_2DRR_Vib_Cas2(:,1))./abs(CH2Cl2_2DRR_Vib_Dir(:,1)),'b-','Linewidth',2);hold on
xlabel('SOLVENT VIBRATIONAL FREQUENCY (cm^{-1})');ylabel('RATIO');
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold');axis square;
legend('FSRS','2DRR');
hold off