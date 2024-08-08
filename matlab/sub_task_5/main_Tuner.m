clear;
% close all;
%
%
%
nmode=1;
wvib(1)=670;
% disp(1)=0.35;
% wvib(2)=1100;
% disp(2)=.35;
% wvib(3)=1500;
% disp(3)=1;
weg=20000;
%
%
%
w_aps=0:100:1000;
disps=0:0.25:1;
for iw=1:length(w_aps)
    for id=1:length(disps)
        w_ap=w_aps(iw)+weg;
        disp=disps(id);

        parameters_material.disp=disp;
        parameters_material.gamma_eg=1000;
        parameters_material.gamma_vib=10;
        parameters_material.weg=20000;
        parameters_material.wvib=1370;
        
        % Consolidate laser parameters into structure variable parameters_laser
        parameters_laser.LAMBDA_ap=210;
        parameters_laser.LAMBDA_rp=40;
        parameters_laser.w_ap=w_ap;
        parameters_laser.w_rp=w_ap;

        
%         [base,wviball]=basis(nmode,wvib);
%         [fcall] = fcinfo(base,disp,nmode);
        [base,wviball]=basis(nmode,wvib);
        [fcall] = fcinfo(base,disp,nmode);

%fcall=abs(fcall);
%         [direct,seq,e3,tau,wrs] = response2_Tuner(wviball,fcall,wvib,w_ap);
%         [direct,seq,e3,tau,wrs] = FSRS_TC_Tuning(wviball,fcall,...
%             3,parameters_material,parameters_laser);
        [d,s,e,~,~] = response2_Tuner(wviball,fcall,wvib,w_ap);
        [dtc,stc,etc,~,~] = FSRS_TC_Tuning(wviball,fcall,...
            4,parameters_material,parameters_laser);
        %
        %
%         ratio_Tuner_AM(id,iw)=3e13*abs(seq)./abs(direct)/3e10/2;
        rAM(id,iw)=3e13*abs(s(100))./abs(d(100))/3e10/2;
        rTC(id,iw)=3e13*abs(stc(100))./abs(dtc(100))/3e10/2;

    end
end

if false
         ratio_Tuner_AM(isnan(ratio_Tuner_AM))=0;
 figure;contour(w_aps,disps,abs(ratio_Tuner_AM),15);
%     axis square;
    colormap jet;colorbar;
    set(gca,'linewidth',2,'fontsize',16);caxis([0 0.1]);
    
     for ii=1:5:100
     ratio_Tuner_AM(isnan(ratio_Tuner_AM))=0;
    contour(w_aps,disps,abs(ratio_Tuner_AM(:,:,ii)),15);
%     axis square;
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);title(num2str(wrs(ii)));
    pause
    end
        ratio_Tuner_AM(isnan(ratio_Tuner_AM))=0;
    contour(w_aps,disps,abs(trapz(ratio_Tuner_AM,3)),15);
%     axis square;
    colormap jet;colorbar

end