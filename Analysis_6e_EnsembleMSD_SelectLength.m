clearvars; 

namestr = '161101_3D3_P11_2d'; 
fitlength = 30; %number of points to fit for MSD, 
minlength = 100; %min MSD to be included
maxlength = 3000; %max MSD to be included 

pixelsize = 71; %nm from microscope
sperframe = 0.1; %s per frame usually 0.1s (10 Hz)
xlimit = [0 10];

load(['Working/' namestr '_msd_new.mat']); 
load(['Working/' namestr '_sizeinter.mat']); 

trajgran = ids(2,:) == 1; 
trajscrum = ids(2,:) == 2;
trajlength = (ids(3,:) >= minlength) & (ids(3,:) <= maxlength); %select the proper length range 

idxgran = find(trajgran & trajlength);
idxscrum = find(trajscrum & trajlength); 


msdgran = msd(:,idxgran);
msdscrum = msd(:,idxscrum);

cpixel = pixelsize^2*(1/1000)^2; %pixel conversion to um^2/s

avgmsdgran = sum(msdgran,2) ./ sum(msdgran>0,2) .* cpixel;
avgmsdscrum = sum(msdscrum,2) ./ sum(msdscrum>0,2) .* cpixel;


t = 1:size(msd,1)/2;
t = sperframe.*t; %convert from frames to s 

%fit the MSD
logt = log10(t(1:fitlength));
opt = optimset('Display','none','TolFun',10^-8,'TolX',10^-8);
lfit = @(x,t)x(1).*t+x(2);
logmsd = log10(avgmsdgran(1:fitlength));
[resultgran,~,res,~,~,~,J] = lsqcurvefit(lfit,[1,1],logt,logmsd',[],[],opt);
ci = nlparci(resultgran,res,'jacobian',J);
plusminusgran(1,:) = ci(1,:) - resultgran(1);
plusminusgran(2,:) = ci(2,:) - resultgran(2);

result = resultgran;
temp(1) = result(1);
temp(2) = ci(1,1);
temp(3) = ci(1,2);
temp(4) = result(2);
temp(5) = ci(2,1);
temp(6) = ci(2,2);
grancol = temp;

logmsd = log10(avgmsdscrum(1:fitlength));
[resultscrum,~,res,~,~,~,J] = lsqcurvefit(lfit,[1,1],logt,logmsd',[],[],opt);
ci = nlparci(resultscrum,res,'jacobian',J);
plusminusscrum(1,:) = ci(1,:) - resultscrum(1);
plusminusscrum(2,:) = ci(2,:) - resultscrum(2);

result = resultscrum;
temp(1) = result(1);
temp(2) = ci(1,1);
temp(3) = ci(1,2);
temp(4) = result(2);
temp(5) = ci(2,1);
temp(6) = ci(2,2);
scrumcol = temp;

fitgran = 10^resultgran(2).*t.^resultgran(1);
fitscrum = 10^resultscrum(2).*t.^resultscrum(1);

figure(1)
subplot(2,1,1)
plot(t,avgmsdgran(1:length(t)),'o',t,fitgran,'-')
title(['Granules: ' num2str(minlength) '-' num2str(maxlength) ...
    ' \alpha = ' num2str(resultgran(1),'%1.2f') ' +/- ' num2str(plusminusgran(1,2),'%1.2f')])
xlim(xlimit)
xlabel('delta (s)')
ylabel('MSD (um^2)')
set(gca,'FontSize',16)
subplot(2,1,2)
plot(t,avgmsdscrum(1:length(t)),'o',t,fitscrum,'-')
title(['Scrum: ' num2str(minlength) '-' num2str(maxlength)...
    ' \alpha = ' num2str(resultscrum(1),'%1.2f') ' +/- ' num2str(plusminusscrum(1,2),'%1.2f')])
xlim(xlimit)
xlabel('delta (s)')
ylabel('MSD (um^2)')
set(gca,'FontSize',16)

dataT = table(t,avgmsdgran(1:length(t))',avgmsdscrum(1:length(t))',fitgran,fitscrum, ...
    'VariableNames',{'time_s','granule_msd','scrum_msd','granule_fit','scrum_fit'});

fitT = table(grancol,scrumcol, ...
    'VariableNames',{'granule_fit','scrum_fit'});

writetable(dataT,['Output/' namestr '_length_' num2str(minlength) '_' num2str(maxlength) ...
    '_fit' num2str(fitlength) '_ensemblefit.csv']);

writetable(fitT,['Output/' namestr '_length_' num2str(minlength) '_' num2str(maxlength) ...
    '_fit' num2str(fitlength) '_ensemblefitparams.csv']);

save(['Working/' namestr '_length_' num2str(minlength) '_' num2str(maxlength) '_ensemblefit.mat'], ...
    't','avgmsdgran','fitgran','resultgran','plusminusgran', ...
    'avgmsdscrum','fitscrum','resultscrum','plusminusscrum', ...
    'fitlength','minlength','maxlength');

