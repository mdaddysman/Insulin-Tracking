clearvars;

namestr = '161101_3D3_P11_2d';
nptsfit = 20;
ms = 16; %markersize 
lw = 2; %linewidth

load([namestr '_eamsd.mat']); 

%fits for figure one
fitt = eatrmsd(2:nptsfit+1,1);
opt = optimset('Display','none','TolFun',10^-8,'TolX',10^-8);
lfit = @(x,t)x(2).*t.^x(1); 

result = zeros(4,3);

for m=1:4
    fitmsd = eatrmsd(2:nptsfit+1,m+1);
    [r,~,res,~,~,~,J] = lsqcurvefit(lfit,[-0.3,1],fitt,fitmsd,[],[],opt);
    result(m,1:2) = r;
    ci = nlparci(r,res,'jacobian',J);
    result(m,3) = ci(1,2) - result(m,1);
end

figure(1)
ph(1) = loglog(eatrmsd(:,1),eatrmsd(:,2),'.','MarkerSize',ms);
rT(1) = corr2(log10(fitt),log10(eatrmsd(2:nptsfit+1,2)));
hold on
for m=3:5
   ph(m-1) = loglog(eatrmsd(:,1),eatrmsd(:,m),'.','MarkerSize',ms);
   rT(m-1) = corr2(log10(fitt),log10(eatrmsd(2:nptsfit+1,m)));
end
for m=1:4
    cph = loglog(eatrmsd(:,1),lfit(result(m,1:3),eatrmsd(:,1)),'-','LineWidth',lw);
    set(cph,'Color',get(ph(m),'Color'));
end
hold off
xlim([0 10])
xlabel('T (s)')
ylabel('MSD (um^2)')
legend({['s=1 T^{' num2str(result(1,1),'%1.3f') ' +/- ' num2str(result(1,3),'%1.3f') '} R = ' num2str(rT(1),'%1.3f')], ...
    ['s=2 T^{' num2str(result(2,1),'%1.3f') ' +/- ' num2str(result(2,3),'%1.3f') '} R = ' num2str(rT(2),'%1.3f')], ...
    ['s=3 T^{' num2str(result(3,1),'%1.3f') ' +/- ' num2str(result(3,3),'%1.3f') '} R = ' num2str(rT(3),'%1.3f')], ...
    ['s=4 T^{' num2str(result(4,1),'%1.3f') ' +/- ' num2str(result(4,3),'%1.3f') '} R = ' num2str(rT(4),'%1.3f')]}, ...
    'Location','Southwest','FontSize',18)
set(gca,'FontSize',16)


figure(2)
semilogy(eamsd(:,1),eamsd(:,2),'.')
r(1) = corr2(eamsd(:,1)+eamsd(2,1),log10(eamsd(:,2)));
hold on
for m=3:5
    semilogy (eamsd(:,1),eamsd(:,m),'.')
    r(m-1) = corr2(eamsd(:,1)+eamsd(2,1),log10(eamsd(:,m)));
end
hold off
%xlim([0 10])
xlabel('time (s)')
ylabel('MSD (um^2)')
legend(['s=1 log-lin R = ' num2str(r(1),'%1.3f')], ...
    ['s=2 log-lin R = ' num2str(r(2),'%1.3f')], ...
    ['s=3 log-lin R = ' num2str(r(3),'%1.3f')], ...
    ['s=4 log-lin R = ' num2str(r(4),'%1.3f')])
set(gca,'FontSize',16)

figure(3)
plot(eamsd(:,1),eamsd(:,2)-mean(eamsd(:,2)),'.')
hold on
for m=3:5
    plot(eamsd(:,1),eamsd(:,m)-mean(eamsd(:,m)),'.')
end
hold off
%xlim([0 10])
xlabel('time (s)')
ylabel('MSD detrend (um^2)')
legend(['s=1 log-log R: ' num2str(r(1))], ...
    ['s=2 log-log R: ' num2str(r(2))], ...
    ['s=3 log-log R: ' num2str(r(3))], ...
    ['s=4 log-log R: ' num2str(r(4))])
set(gca,'FontSize',16)

autocorr = zeros(size(eamsd,1)*2-1,4);
crosscorr = zeros(size(eamsd,1)*2-1,4);

for m=1:4
    autocorr(:,m) = xcorr(eamsd(:,m+1)-mean(eamsd(:,m+1)),'coeff');
    crosscorr(:,m) = xcorr(eamsd(:,m+1)-mean(eamsd(:,m+1)),eamsd(:,2)-mean(eamsd(:,2)),'coeff');
end

t = -eamsd(end,1):eamsd(2,1):eamsd(end,1);

figure(4)
plot(t,autocorr(:,1),'-o')
hold on
for m=2:4
    plot(t,autocorr(:,m),'-o')
end
hold off
%xlim([0 10])
xlabel('time offset (s)')
ylabel('autocorrleation')
legend('s=1','s=2','s=3','s=4')
set(gca,'FontSize',16)

figure(5)
for m=1:4
    subplot(4,1,m)
    tph = plot(t,autocorr(:,m),'-o');
    ylabel('autocorrleation')
    title(['s = ' num2str(m)])
    set(gca,'FontSize',16)
    ylim([-0.05 1])
    set(tph,'Color',get(ph(m),'Color'));
end
%xlim([0 10])
xlabel('time offset (s)')
tightfig;



figure(6)
plot(t,crosscorr(:,1),'-o')
hold on
for m=2:4
    plot(t,crosscorr(:,m),'-o')
end
hold off
%xlim([0 10])
xlabel('time offset (s)')
ylabel('crosscorrleation')
legend('s=1xs=1','s=1xs=2','s=1xs=3','s=1xs=4')
set(gca,'FontSize',16)