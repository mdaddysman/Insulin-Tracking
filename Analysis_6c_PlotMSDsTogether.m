clearvars;

namestr = {'161101_3D3_P11_2b_bottom'; ...
    '161101_3D3_P11_2d'; ...
    '161101_3D3_P11_2c_top'};

xlimit = [0 30];
ms = 26;
lw = 3;

alphaboth = zeros(length(namestr),2);
alphaclass = zeros(length(namestr),2);

figure('Name','Granules of longer than fit length')
for m=1:length(namestr)
    load([namestr{m} '_class' num2str(1) '_ensemblefit.mat'])
    avgmsdboth1 = avgmsdboth;
    fitboth1 = fitboth;
    t1 = t;
    alphaboth(m,1) = resultboth(1);
    load([namestr{m} '_class' num2str(2) '_ensemblefit.mat'])
    alphaboth(m,2) = resultboth(1);
    subplot(3,1,m)
    ph =plot(t1,avgmsdboth1(1:length(t1)),'.',t1,fitboth1,'-',t,avgmsdboth(1:length(t)),'.',t,fitboth,'-', ...
        'MarkerSize',ms,'LineWidth',lw);
    xlim(xlimit)
    maxidx = find(t1 == xlimit(2));
    ymin = min([fitboth1(1) fitboth]);
    ymax = max([fitboth1(maxidx) fitboth(maxidx)]);
    ylim([ymin ymax])
    set(gca,'FontSize',16,'XScale','log','YScale','log');
    if(m ~= 3)
        set(gca,'XLabel',[]);
    end
    ph(1).Color = [0.216 0.494 0.722];
    ph(2).Color = [0.216 0.494 0.722];
    ph(3).Color = [0.894 0.102 0.11];
    ph(4).Color = [0.894 0.102 0.11];
end

figure('Name','Granules of all trajectory lengths')
for m=1:length(namestr)
    load([namestr{m} '_class' num2str(1) '_ensemblefit.mat'])
    avgmsdclass1 = avgmsdclass;
    fitclass1 = fitclass;
    t1 = t;
    alphaclass(m,1) = resultclass(1);
    load([namestr{m} '_class' num2str(2) '_ensemblefit.mat'])
    alphaclass(m,2) = resultclass(1);
    subplot(3,1,m)
    ph = plot(t1,avgmsdclass1(1:length(t1)),'.',t1,fitclass1,'-',t,avgmsdclass(1:length(t)),'.',t,fitclass,'-', ...
        'MarkerSize',ms,'LineWidth',lw);
    xlim(xlimit)
    maxidx = find(t1 == xlimit(2));
    ymin = min([fitclass1(1) fitclass]);
    ymax = max([fitclass1(maxidx) fitclass(maxidx)]);
    ylim([ymin ymax])
    set(gca,'FontSize',16,'XScale','log','YScale','log');
    if(m ~= 3)
        set(gca,'XLabel',[]);
    end
    ph(1).Color = [0.216 0.494 0.722];
    ph(2).Color = [0.216 0.494 0.722];
    ph(3).Color = [0.894 0.102 0.11];
    ph(4).Color = [0.894 0.102 0.11];
end

figure('Name','Granules of longer than fit length, single plot log')
for m=1:length(namestr)
    load([namestr{m} '_class' num2str(1) '_ensemblefit.mat'])
    avgmsdboth1 = avgmsdboth;
    fitboth1 = fitboth;
    t1 = t;
    alphaboth(m,1) = resultboth(1);
    pmboth(m,1:2) = plusminusboth(1,1:2);
    load([namestr{m} '_class' num2str(2) '_ensemblefit.mat'])
    alphaboth(m,2) = resultboth(1);
    pmboth(m,3:4) = plusminusboth(1,1:2);
    subplot(2,2,m)
    ph = plot(t1,avgmsdboth1(1:length(t1)),'.',t1,fitboth1,'-',t,avgmsdboth(1:length(t)),'.',t,fitboth,'-', ...
        'MarkerSize',ms,'LineWidth',lw);
    xlim(xlimit)
    maxidx = find(t1 == xlimit(2));
    ymin = min([fitboth1(1) fitboth]);
    ymax = max([fitboth1(maxidx) fitboth(maxidx)]);
    ylim([ymin ymax])
    set(gca,'FontSize',16,'XScale','log','YScale','log','LineWidth',2,'TickLength',3.*[0.01 0.025]);
    ph(1).Color = [0.216 0.494 0.722];
    ph(2).Color = [0.216 0.494 0.722];
    ph(3).Color = [0.894 0.102 0.11];
    ph(4).Color = [0.894 0.102 0.11];
    ylabel('MSD (um^2)')
    xlabel('delta (s)')
    title(num2str(m),'FontSize',8)
end
tightfig;