[~,midx] = max(ids(3,:));
%midx = find(ids(3,:) > 1000);
midx = midx(1);
fidx = find(ID == ids(1,midx));

ltime = log10(0.1.*msd(1:ids(3,midx)/2,1));

fitline = alpha(fidx).*ltime+D(fidx);

figure(1)
plot(ltime,log10(msd(1:ids(3,midx)/2,midx)),'o',ltime,fitline,'LineWidth',2)
xlabel('log10 time')
ylabel('log10 MSD')
set(gca,'FontSize',16)