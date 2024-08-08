clear;
% close all;
%
%
%
nmode=1;
wvib(1)=1370;
disp(1)=0.35;
% wvib(2)=1100;
% disp(2)=.35;
% wvib(3)=1500;
% disp(3)=1;
%
%
%
wvibs=000:100:2000;
disps=0.1:0.1:1;
for iw=1:length(wvibs)
    for id=1:length(disps)
        wvib=wvibs(iw);
        disp=disps(id);

        [base,wviball]=basis(nmode,wvib);
        [fcall] = fcinfo(base,disp,nmode);
        %fcall=abs(fcall);
        [direct,seq,e3,tau,wrs] = response2_looper(wviball,fcall,wvib);
        %
        %
        ratio_AM(id,iw,:)=3e13*abs(seq(100))./abs(direct(100))/3e10/2;
    end
end

if false
     ratio_AM(isnan(ratio_AM))=0;
    figure;contour(disps,wvibs,ratio_AM,35);%     axis square;
    colormap jet;colorbar
    set(gca,'linewidth',2,'fontsize',16);
    caxis([0.0 0.1]);

end