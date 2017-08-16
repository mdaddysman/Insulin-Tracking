namestr = '161101_3D3_P11_2d';
fontsize = 16;

dataset = 's'; % s or l for small or large 

if strcmp(dataset,'s') == 1
    load([namestr '_SPIFF.mat'],'smalldata');
    clist = smalldata(:,2:3);
    load([namestr '_filtered.mat'],'smalldata');
    data = smalldata(:,2:3);
else
    load([namestr '_SPIFF.mat'],'largedata');
    clist = largedata(:,2:3);
    load([namestr '_filtered.mat'],'largedata');
    data = largedata(:,2:3);
end
    
mp_totx = data(:,1)-floor(data(:,1));
mp_toty = data(:,2)-floor(data(:,2));
corr_totx = clist(:,1) - floor(clist(:,1));
corr_toty = clist(:,2) - floor(clist(:,2));
correction = data - (clist-0.5);
corr2 = [mp_totx mp_toty] - [corr_totx corr_toty];

dat = [mp_totx mp_toty];
h = figure;
n = hist3(dat,[30 30]); 
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
close(h)

figure(1)
subplot(1,2,1)
plot(mp_totx,mp_toty,'.'); 
xlabel('Meta-Pixel location','Fontsize',fontsize); 
ylabel('Meta-Pixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
set(gca,'FontSize',fontsize);
axis square
subplot(1,2,2)
pcolor(xb,yb,n1);
xlabel('Meta-Pixel location','Fontsize',fontsize); 
ylabel('Meta-Pixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
colorbar;
set(gca,'FontSize',fontsize);
set(gcf,'Position',[17,574,1501,665]);


dat2 = [corr_totx corr_toty];
h = figure;
n = hist3(dat2,[30 30]); 
n2 = n';
n2(size(n,1) + 1, size(n,2) + 1) = 0;
xb2 = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb2 = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
close(h)

figure(2)
subplot(1,2,1)
plot(corr_totx,corr_toty,'.'); 
xlabel('Meta-Pixel location','Fontsize',fontsize); 
ylabel('Meta-Pixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
set(gca,'FontSize',fontsize);
subplot(1,2,2)
pcolor(xb2,yb2,n2);
xlabel('Meta-Pixel location','Fontsize',fontsize); 
ylabel('Meta-Pixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
colorbar;
set(gca,'FontSize',fontsize);
set(gcf,'Position',[17,574,1501,665]);


figure(3)

plot(correction(:,1),correction(:,2),'.'); axis equal


figure(4)
plot(corr2(:,1),corr2(:,2),'.')

xshift = 0.1;
yshift = 0.060;
lw = 3;
tl = 0.02;
figure(5)
subplot(2,2,1)
plot(mp_totx,mp_toty,'.'); 
xlabel('subpixel location','Fontsize',fontsize); 
ylabel('subpixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
set(gca,'FontSize',fontsize, 'XTick', [0 0.25 0.50 0.75 1], 'YTick', [0 0.25 0.50 0.75 1],'LineWidth',lw, ...
    'Layer','top','TickLength',[tl,tl]);
subplot(2,2,2)
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
xlabel('subpixel location','Fontsize',fontsize); 
ylabel('subpixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
colorbar;
set(gca,'FontSize',fontsize, 'XTick', [0 0.25 0.50 0.75 1], 'YTick', [0 0.25 0.50 0.75 1],'LineWidth',lw, ...
    'Layer','top','TickLength',[tl,tl]);
p = get(gca,'Position');
set(gca,'Position',[p(1) - xshift, p(2:4)]);
subplot(2,2,3)
plot(corr_totx,corr_toty,'.'); 
xlabel('subpixel location','Fontsize',fontsize); 
ylabel('subpixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
set(gca,'FontSize',fontsize, 'XTick', [0 0.25 0.50 0.75 1], 'YTick', [0 0.25 0.50 0.75 1],'LineWidth',lw, ...
    'Layer','top','TickLength',[tl,tl]);
p = get(gca,'Position');
set(gca,'Position',[p(1), p(2) + yshift, p(3:4)]);
subplot(2,2,4)
h = pcolor(xb2,yb2,n2);
set(h, 'EdgeColor', 'none');
xlabel('subpixel location','Fontsize',fontsize); 
ylabel('subpixel location','Fontsize',fontsize); 
xlim([0 1]); ylim([0 1]);
axis square
colorbar;
set(gca,'FontSize',fontsize, 'XTick', [0 0.25 0.50 0.75 1], 'YTick', [0 0.25 0.50 0.75 1],'LineWidth',lw, ...
    'Layer','top','TickLength',[tl,tl]);
p = get(gca,'Position');
set(gca,'Position',[p(1) - xshift, p(2) + yshift, p(3:4)]);