clearvars;

namestr = '161101_3D3_P11_2d';
lw = 3;
ms = 20; 
limity = [10^-4 10^-1];
limitx = [-200 200];
pixelsize = 71; %nm 

load([namestr '_sizeinter.mat']);

stracked = tracked(tracked(:,5) == 1,1:5);
ltracked = tracked(tracked(:,5) == 2,1:5);

stracked = sortrows(stracked,[4 3]);
ltracked = sortrows(ltracked,[4 3]);

sstep = stracked(2:end,:) - stracked(1:end-1,:);
sstep = sstep(sstep(:,4) == 0,:);

lstep = ltracked(2:end,:) - ltracked(1:end-1,:);
lstep = lstep(lstep(:,4) == 0,:);

sstep2 = sstep(:,1:2).*pixelsize;
lstep2 = lstep(:,1:2).*pixelsize;


sxhist = histogrampts(sstep2(:,1));
syhist = histogrampts(sstep2(:,2));

lxhist = histogrampts(lstep2(:,1));
lyhist = histogrampts(lstep2(:,2));

resultsx = lsqcurvefit('gauss1D_noC',[1 0 100],sxhist(:,1),sxhist(:,2));
resultsy = lsqcurvefit('gauss1D_noC',[1 0 100],syhist(:,1),syhist(:,2));

resultlx = lsqcurvefit('gauss1D_noC',[1 0 100],lxhist(:,1),lxhist(:,2));
resultly = lsqcurvefit('gauss1D_noC',[1 0 100],lyhist(:,1),lyhist(:,2));

t = -1000:1000;

titlestr = strrep(namestr,'_',' ');

figure(1)
subplot(1,2,1)
semilogy(sxhist(:,1),sxhist(:,2),'o',t,gauss1D_noC(resultsx,t),'-')
ylim(limity);
xlim(limitx);
set(gca,'FontSize',16)
subplot(1,2,2)
semilogy(syhist(:,1),syhist(:,2),'o',t,gauss1D_noC(resultsy,t),'-')
ylim(limity);
xlim(limitx);
set(gca,'FontSize',16)

figure(2)
subplot(1,2,1)
semilogy(lxhist(:,1),lxhist(:,2),'o',t,gauss1D_noC(resultlx,t),'-')
ylim(limity);
xlim(limitx);
set(gca,'FontSize',16)
subplot(1,2,2)
semilogy(lyhist(:,1),lyhist(:,2),'o',t,gauss1D_noC(resultly,t),'-')
ylim(limity);
xlim(limitx);
set(gca,'FontSize',16)


figure(3)
semilogy(sxhist(:,1),sxhist(:,2),'.',t,gauss1D_noC(resultsx,t),'-','LineWidth',lw,'MarkerSize',ms)
ylim(limity);
xlim(limitx);
set(gca,'FontSize',16)
title([titlestr ' Small X Step'])
xlabel('x step (nm)')
ylabel('pdf')
legend('data','gauss fit')

figure(4)
semilogy(lxhist(:,1),lxhist(:,2),'.',t,gauss1D_noC(resultlx,t),'-','LineWidth',lw,'MarkerSize',ms)
ylim(limity);
xlim(limitx);
set(gca,'FontSize',16)
title([titlestr ' Large X Step'])
xlabel('x step (nm)')
ylabel('pdf')
legend('data','gauss fit')


