namestr = 'cMovie1';
pixelsize = 71; %in nm
fontsize = 16;
pixelsperbin = 10;
filtFWHM = 1.5; %in pixels

load([namestr '_SPIFF.mat']);

%Plot all the locations
figure(1)
set(gcf,'Name','Positions');
plot(smalldata(:,2),smalldata(:,3),'.',largedata(:,2),largedata(:,3),'.');
legend('Granules','Scrums');
set(gca,'FontSize',fontsize)

%Generate a histogram of the locations
%granules
xlimits = [min(smalldata(:,2)) max(smalldata(:,2))];
ylimits = [min(smalldata(:,3)) max(smalldata(:,3))];
xrange = xlimits(2) - xlimits(1);
yrange = ylimits(2) - ylimits(1);
nbins = round([xrange yrange]./pixelsperbin);
h = figure;
sn = hist3(smalldata(:,2:3),nbins); 
xsb = linspace(xlimits(1),xlimits(2),size(sn,1));
ysb = linspace(ylimits(1),ylimits(2),size(sn,2));
close(h)

%scrums
xlimits = [min(largedata(:,2)) max(largedata(:,2))];
ylimits = [min(largedata(:,3)) max(largedata(:,3))];
xrange = xlimits(2) - xlimits(1);
yrange = ylimits(2) - ylimits(1);
nbins = round([xrange yrange]./pixelsperbin);
h = figure;
ln = hist3(largedata(:,2:3),nbins); 
xlb = linspace(xlimits(1),xlimits(2),size(ln,1));
ylb = linspace(ylimits(1),ylimits(2),size(ln,2));
close(h)

figure(2)
pcolor(xsb,ysb,sn');
set(gcf,'Name','Granule Histogram');
colormap inferno
colorbar;
set(gca,'FontSize',fontsize,'Color','k');
xlim([0 ceil(xsb(end)/100)*100]);
ylim([0 ceil(ysb(end)/100)*100]);

figure(3)
pcolor(xlb,ylb,ln');
set(gcf,'Name','Scrum Histogram');
colormap inferno
colorbar;
set(gca,'FontSize',fontsize,'Color','k');
xlim([0 ceil(xsb(end)/100)*100]);
ylim([0 ceil(ysb(end)/100)*100]);


%make a single pixel histogram for the granules
%smooth by a FWHM gaussian from start 
%normalize and plot it 
xlimits = [min(smalldata(:,2)) max(smalldata(:,2))];
ylimits = [min(smalldata(:,3)) max(smalldata(:,3))];
xrange = xlimits(2) - xlimits(1);
yrange = ylimits(2) - ylimits(1);
nbins = ceil([xrange yrange]);
h = figure;
nsb = hist3(smalldata(:,2:3),nbins); 
xdb = linspace(xlimits(1),xlimits(2),size(nsb,1));
ydb = linspace(ylimits(1),ylimits(2),size(nsb,2));
close(h)

sigfilt = filtFWHM/(2*sqrt(2*log(2)));%convert FWHM to sigma
nsb_flit = imgaussfilt(nsb,sigfilt);
nsb_flit_norm = nsb_flit./max(max(nsb_flit));

figure(4)
imshow(flipud(255.*nsb_flit_norm'),inferno);
set(gcf,'Name','Smoothed and Normalized Histogram');
colorbar;
set(gca,'FontSize',fontsize);

nsb_thres = nsb_flit_norm>0.01;

figure(5)
imshow(flipud(nsb_thres'));
set(gcf,'Name','Thresholded');