

% Load data from mat file
load 2dRR3
set(0,'units','normalized');
dw=1/nt;c=3e-5;
maxr=max(max(max(abs(ratio))));
% Axes positions
gcpos=[0.05 0.605 0.3 0.3; 0.05 0.14 0.3 0.3];
gdpos=[0.30 0.605 0.3 0.3; 0.3 0.14 0.3 0.3];
grpos=[0.55 0.605 0.3 0.3; 0.55 0.14 0.3 0.3];
% gcpos=[0.05 0.59 0.33 0.33; 0.05 0.09 0.33 0.33];
% gdpos=[0.33 0.59 0.33 0.33; 0.33 0.09 0.33 0.33];
% grpos=[0.61 0.59 0.33 0.33; 0.61 0.09 0.33 0.33];
% Color bar positions
ccpos=[0.115 0.92 0.172 0.03;0.115 0.455 0.172 0.03];
cdpos=[0.365 0.92 0.172 0.03;0.365 0.455 0.172 0.03];
crpos=[0.615 0.92 0.172 0.03;0.615 0.455 0.172 0.03];
% ccpos=[0.126 0.925 0.178 0.03;0.126 0.425 0.178 0.03];
% cdpos=[0.406 0.925 0.178 0.03;0.406 0.425 0.178 0.03];
% crpos=[0.686 0.925 0.178 0.03;0.686 0.425 0.178 0.03];
% Text
tstr={'a)';'b)';'c)';'d)';'e)';'f)'};
% Text positions
tpos=[0.26 0.84 0.06 0.06;0.51 0.84 0.06 0.06;0.76 0.84 0.06 0.06;...
    0.26 0.37 0.06 0.06;0.51 0.37 0.06 0.06;0.76 0.37 0.06 0.06];
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 1 1]);%[0 4/75 1 0.85]);
for iv=1:2
    maxs=max(max(max(prefactor(iv)/3e10*abs(cascade2d(:,:,iv)))),...
        max(max(abs(direct2d(:,:,iv)))));
%     dt=dts(iv);
    w=(-1/2:dw:(1/2-dw))/dt/c;
    gc(iv)=subplot(2,3,1+3*(iv-1));contour(w,w,prefactor(iv)/3e10*abs(cascade2d(:,:,iv))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');xlim([-2000 2000]);
    ylabel('\omega_2/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;cc(iv)=colorbar('Location','northoutside','fontsize',10,'Position',ccpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'Position',gcpos(iv,:));
    axis square;%title('CASCADES')
    tx=annotation(gcf,'textbox',tpos(1+3*(iv-1),:),...
        'String',tstr(1+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);
    
    gd(iv)=subplot(2,3,2+3*(iv-1));contour(w,w,abs(direct2d(:,:,iv))/maxs,50);
    xlabel('\omega_1/2\pic (cm^{-1})');xlim([-2000 2000]);
    ylabel('\omega_2/2\pic (cm^{-1})');ylim([-2000 2000]);
    colormap jet;cd(iv)=colorbar('Location','northoutside','fontsize',10,'Position',cdpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'Position',gdpos(iv,:));axis square;%title('Direct')
    tx=annotation(gcf,'textbox',tpos(2+3*(iv-1),:),...
        'String',tstr(2+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    gr(iv)=subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})');
    ylabel({'Mode';'Displacement'});caxis([0 maxr])
    colormap jet;cr(iv)=colorbar('Location','northoutside','fontsize',10,'Position',crpos(iv,:));
    set(gca,'linewidth',2,'fontsize',10,'Position',grpos(iv,:));axis square;%title('E^{3}_{CAS}:E^{5}_{DIRECT}')
    tx=annotation(gcf,'textbox',tpos(3+3*(iv-1),:),...
        'String',tstr(3+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);
end