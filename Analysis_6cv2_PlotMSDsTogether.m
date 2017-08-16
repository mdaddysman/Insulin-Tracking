clearvars;

namestr = '161101_3D3_P11_2d';

minlengths = [10 10 100 1000];
maxlengths = [3000 100 1000 3000];

xlimit = [0 30];
ms = 26;
lw = 3;

ymin = 1*10^-4;
ymax = 2*10^1;

figure('Name',strrep(namestr,'_',' '))
for m=1:length(minlengths)
    load([namestr '_length_' num2str(minlengths(m)) '_' num2str(maxlengths(m)) '_ensemblefit.mat'])
    subplot(4,1,m)
     
    ph = plot(t,avgmsdgran(1:length(t)),'.',t,fitgran,'-',t,avgmsdscrum(1:length(t)),'.',t,fitscrum,'-', ...
        'MarkerSize',ms,'LineWidth',lw);
    xlim(xlimit)
    maxidx = find(t == xlimit(2));
    %ymin = min([fitgran(1) fitscrum(1)]);
    %ymax = max([fitgran(maxidx) fitscrum(maxidx)]);
    ylim([ymin ymax])
    set(gca,'FontSize',16,'XScale','log','YScale','log','LineWidth',2,'TickLength',3.*[0.01 0.025]);
    ph(1).Color = [0.216 0.494 0.722];
    ph(2).Color = [0.216 0.494 0.722];
    ph(3).Color = [0.894 0.102 0.11];
    ph(4).Color = [0.894 0.102 0.11];
    ylabel('MSD (\mum^2)')
    if m == length(minlengths)
        xlabel('\Delta (s)')
    end
    title(['Length: ' num2str(minlengths(m)) '-' num2str(maxlengths(m))])
    
    text(t(2),ymax-0.85*ymax,['\alpha = ' num2str(resultgran(1),'%1.2f') ' +/- ' num2str(plusminusgran(1,2),'%1.2f')], ...
        'Color', [0.216 0.494 0.722],'FontSize',16);
    text(t(2),ymax-0.97*ymax,['\alpha = ' num2str(resultscrum(1),'%1.2f') ' +/- ' num2str(plusminusscrum(1,2),'%1.2f')], ...
        'Color', [0.894 0.102 0.11],'FontSize',16);
    patch('XData',[t(1) t(1) t(fitlength) t(fitlength)],'YData', [ymin ymax ymax ymin],...
        'FaceColor', 0.7.*ones(3,1),'EdgeColor','none','FaceAlpha',0.5);
    set(gca,'children',flipud(get(gca,'children')))
end
tightfig;
set(gcf,'Position',[138 331 440 757]);