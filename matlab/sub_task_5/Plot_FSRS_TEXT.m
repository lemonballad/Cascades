% Load data from mat file
load FSRS2
set(0,'units','normalized');
maxr=max(max(max(abs(ratio))));
% Axes positions
gcpos=[0.05 0.605 0.3 0.3; 0.05 0.14 0.3 0.3];
gdpos=[0.30 0.605 0.3 0.3; 0.3 0.14 0.3 0.3];
grpos=[0.55 0.605 0.3 0.3; 0.55 0.14 0.3 0.3];
% Color bar positions
ccpos=[0.115 0.92 0.172 0.03;0.115 0.455 0.172 0.03];
cdpos=[0.365 0.92 0.172 0.03;0.365 0.455 0.172 0.03];
crpos=[0.615 0.92 0.172 0.03;0.615 0.455 0.172 0.03];
% Text
tstr={'a)';'b)';'c)';'d)';'e)';'f)'};
% Text positions
tpos=[0.26 0.84 0.06 0.06;0.51 0.84 0.06 0.06;0.76 0.84 0.06 0.06;...
    0.26 0.37 0.06 0.06;0.51 0.37 0.06 0.06;0.76 0.37 0.06 0.06];
gcf=figure;
set(gcf,'Units','Normalized','Position',[0 0 1 1]);%[0 4/75 1 0.85]);
for iv=1:2
    maxs=max(max(max(prefactor(iv)/3e10*abs(cascadeFSRS(:,:,iv)))),...
        max(max(abs(directFSRS(:,:,iv)))));
    
    % Cascades
    subplot(2,3,1+3*(iv-1));contour(t/1e3,w,prefactor(iv)/3e10*abs(cascadeFSRS(:,:,iv))'/maxs,50);
    gc(iv)=get(gca);
    xlabel('\tau_1 (ps)','fontweight','bold');xlim([0 7]);
    ylabel('\omega/2\pic (cm^{-1})','fontweight','bold');ylim([0 2000]);
    %title('CASCADES','fontweight','bold')
    colormap jet;hcb=colorbar('Location','northoutside',...
        'fontweight','bold','fontsize',10,'Position',ccpos(iv,:));
    cc(iv)=get(hcb);
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10,...
        'Position',gcpos(iv,:));axis square;
    tx=annotation(gcf,'textbox',tpos(1+3*(iv-1),:),...
        'String',tstr(1+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    % Direct
    subplot(2,3,2+3*(iv-1));contour(t/1e3,w,abs(directFSRS(:,:,iv))'/maxs,50);
    gd(iv)=get(gca);
    xlabel('\tau_1 (ps)','fontweight','bold');xlim([0 7]);
    ylabel('\omega/2\pic (cm^{-1})','fontweight','bold');ylim([0 2000]);
    %title('Direct','fontweight','bold')
    colormap jet;hcb=colorbar('Location','northoutside',...
        'fontweight','bold','fontsize',10,'Position',cdpos(iv,:));
    cd(iv)=get(hcb);
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10,...
        'Position',gdpos(iv,:));axis square;
    tx=annotation(gcf,'textbox',tpos(2+3*(iv-1),:),...
        'String',tstr(2+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10);

    % Ratio
    subplot(2,3,3+3*(iv-1));contourf(dws,disps,ratio(:,:,iv),50,'edgecolor','none');
    gr(iv)=get(gca);
    xlabel('(\omega_L-\omega_{eg})/2\pic (cm^{-1})','fontweight','bold');
    ylabel({'Mode';'Displacement \itd'},'fontweight','bold');
    %title('E^{3}_{CAS}:E^{5}_{DIRECT}','fontweight','bold')
    colormap jet;hcb=colorbar('Location','northoutside',...
        'fontweight','bold','fontsize',10,'Position',crpos(iv,:));
    cr(iv)=get(hcb);caxis([0 maxr])
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10,...
        'Position',grpos(iv,:));axis square;
    tx=annotation(gcf,'textbox',tpos(3+3*(iv-1),:),...
        'String',tstr(3+3*(iv-1)),'EdgeColor','none',...
        'fontweight','bold','fontsize',10,'Color',[1 1 1]);
end