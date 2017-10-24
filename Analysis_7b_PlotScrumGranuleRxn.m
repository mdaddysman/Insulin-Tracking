clearvars;

namestr = 'cMovie1';
sperframe = 0.1;
nframes = 3000;

radius = 15; %how many pixels must the granule be in to be considered interacting 

load(['Working/' namestr '_radius_' num2str(radius) '_granulescruminter.mat']); 

intscrums = unique(interactions(:,7));
numofintscrums = length(intscrums);

if ishghandle(1)
    close(1)
end

%Enter scrum = end of granule (2)
%Exit scrum = begin of granule (1)

enterevents = interactions(interactions(:,2) == 2,:);
exitevents = interactions(interactions(:,2) == 1,:);

nenterevents = size(enterevents,1);
nexitevents = size(exitevents,1);

enterplot = zeros(nenterevents,2);
exitplot = zeros(nexitevents,2);

enterplot(:,1) = sperframe.*enterevents(:,4);
exitplot(:,1) = sperframe.*exitevents(:,4);

for m=1:nenterevents
    idx = find(enterevents(m,7) == intscrums);
    enterplot(m,2) = idx;
end

for m=1:nexitevents
    idx = find(exitevents(m,7) == intscrums);
    exitplot(m,2) = idx;
end


%used to make fig 1 here 

%figure out when the next enter / exit event occurs after the previous
%event

%let's make a histogram of the length of trajectories first 
figure(2)
histogram(enterevents(:,3).*sperframe,'FaceColor','b','Normalization','pdf','BinMethod','fd')
hold on
histogram(exitevents(:,3).*sperframe,'FaceColor','r','Normalization','pdf','BinMethod','fd')
histogram(potinter(:,3).*sperframe,'FaceColor',0.25.*ones(3,1),'Normalization','pdf','BinMethod','fd')
hold off
legend('Enter','Exit','All') %all includes none interacting trajectories as well
xlabel('length of trajectory (s)')
ylabel('pdf')
xlim([0 60])
title('distribution of trajectory lengths')
set(gca,'FontSize',16);


figure(3)
enterhist = histogram2pts(enterevents(:,3).*sperframe,'Normalization','pdf','BinMethod','fd');
exithist = histogram2pts(exitevents(:,3).*sperframe,'Normalization','pdf','BinMethod','fd');
allhist = histogram2pts(potinter(:,3).*sperframe,'Normalization','pdf','BinMethod','fd');
dph = plot(enterhist(:,1),enterhist(:,2),'-ob',exithist(:,1),exithist(:,2),'-or',allhist(:,1),allhist(:,2),'-o', ...
    'LineWidth',2.5,'MarkerSize',8);
set(dph(3),'Color',0.25.*ones(3,1));
legend('Enter','Exit','All') %all includes none interacting trajectories as well
xlabel('length of trajectory (s)')
ylabel('pdf')
xlim([0 60])
title('distribution of trajectory lengths')
set(gca,'FontSize',16);

%process enter 
ids = unique(enterevents(:,7));
tenternext = [];
for m=1:length(ids)
    temp = squeeze(enterevents(enterevents(:,7) == ids(m),4));
    temp = sort(temp);
    if length(temp) == 1
        continue;
    end
    
    temp2 = temp(2:end) - temp(1:end-1);
    
    tenternext = [tenternext temp2'];
end

%process exit
ids = unique(exitevents(:,7));
texitnext = [];
for m=1:length(ids)
    temp = squeeze(exitevents(exitevents(:,7) == ids(m),4));
    temp = sort(temp);
    if length(temp) == 1
        continue;
    end
    
    temp2 = temp(2:end) - temp(1:end-1);
    
    texitnext = [texitnext temp2'];
end

%process enter -> exit event 
ids_enter = unique(enterevents(:,7));
ids_exit = unique(exitevents(:,7));
ids_both = intersect(ids_enter,ids_exit);
tenterexit = [];
enterexitidstime = [];
for m=1:length(ids_both)
    enters = squeeze(enterevents(enterevents(:,7) == ids_both(m),4));
    enters = sort(enters);
    exits = squeeze(exitevents(exitevents(:,7) == ids_both(m),4));
    exits = sort(exits);
    
    %scroll through the enters and find the next exit and next enter
    for n=1:length(enters)
        temp = exits - enters(n); %find the distance from the current enter to the exits
        temp2 = sort(temp(temp >= 0));
        if isempty(temp2)
            break; %no more exits so break the loop 
        end
        currval = temp2(1);
        exitidx = find(currval == temp,1); 
        if n ~= length(enters) %if it isn't the last enter even check the next one
            if exits(exitidx) > enters(n+1) %check to make sure there isn't another enter event before this exit 
                continue; %catch on the next loop
            end
        end
        %if all this passes we have a enter -> exit event found so record
        %it
        tdiff = exits(exitidx) - enters(n);
        tenterexit = [tenterexit tdiff];
        enterexitidstime = [enterexitidstime; [ids_both(m) enters(n) exits(exitidx)]];
    end
    
end

enter2exitplot = [[enterexitidstime(:,2) enterexitidstime(:,1)]; [enterexitidstime(:,3) enterexitidstime(:,1)]];
enter2exitplot(:,1) = enter2exitplot(:,1).*sperframe;
for m=1:length(enter2exitplot)
    idx = find(enter2exitplot(m,2) == intscrums);
    enter2exitplot(m,2) = idx;
end
figure(1)
ph = plot(enterplot(:,1),enterplot(:,2),'.b',exitplot(:,1),exitplot(:,2),'.r', ...
    enter2exitplot(:,1),enter2exitplot(:,2),'og','MarkerSize',16);
ph(3).LineWidth = 4;
ph(3).MarkerSize = 8;
ph(3).Color = [34, 139, 34]./255;
xlim([0 sperframe*nframes]);
ylim([0 numofintscrums+1]);
xlabel('time (s)')
ylabel('scrum')
legend('Enter','Exit','Location','Northwest')
set(gca,'FontSize',16,'YTick',[])
for m=1:numofintscrums
    temp = scrums(scrums(:,4) == intscrums(m),:);
    patch('XData',[sperframe*temp(1,3) sperframe*temp(1,3) sperframe*temp(end,3) sperframe*temp(end,3)], ...
        'YData', [m-0.25 m+0.25 m+0.25 m-0.25],...
        'FaceColor','black','EdgeColor','none');
end
set(gca,'children',flipud(get(gca,'children')))


tightfig;
set(gcf,'Position',[100 190 975 850]);


tenternext = tenternext.*sperframe;
texitnext = texitnext.*sperframe;
tenterexit = tenterexit.*sperframe;

figure(4)
histogram(tenternext,'FaceColor','b','Normalization','pdf','BinMethod','fd')
hold on
histogram(texitnext,'FaceColor','r','Normalization','pdf','BinMethod','fd')
histogram(tenterexit,'FaceColor','g','Normalization','pdf','BinMethod','fd')
hold off
legend('Enter','Exit','Enter -> Exit') 
xlabel('time until next event (s)')
ylabel('pdf')
title('distribution of time until next event')
set(gca,'FontSize',16);

figure(5)
histogram(tenternext,'FaceColor','b','Normalization','cdf','BinMethod','fd')
hold on
histogram(texitnext,'FaceColor','r','Normalization','cdf','BinMethod','fd')
histogram(tenterexit,'FaceColor','g','Normalization','cdf','BinMethod','fd')
hold off
legend('Enter','Exit','Enter -> Exit') 
xlabel('time until next event (s)')
ylabel('cdf')
title('cdf of time until next event')
set(gca,'FontSize',16);

entertime = histogram2pts(tenternext,'Normalization','pdf','BinMethod','fd');
exittime = histogram2pts(texitnext,'Normalization','pdf','BinMethod','fd');
enter2exittime = histogram2pts(tenterexit,'Normalization','pdf','BinMethod','fd');

opt = optimset('Display','none','TolFun',10^-8,'TolX',10^-8);
efit = @(x,t)x(1).*exp(-x(2).*t);

[r,~,res,~,~,~,J] = lsqcurvefit(efit,[1,1],entertime(:,1),entertime(:,2),[],[],opt);
result(1,1:2) = r;
ci = nlparci(r,res,'jacobian',J);
result(1,3) = ci(2,2) - r(2);

[r,~,res,~,~,~,J] = lsqcurvefit(efit,[1,1],exittime(:,1),exittime(:,2),[],[],opt);
result(2,1:2) = r;
ci = nlparci(r,res,'jacobian',J);
result(2,3) = ci(2,2) - r(2);

[r,~,res,~,~,~,J] = lsqcurvefit(efit,[1,1],enter2exittime(:,1),enter2exittime(:,2),[],[],opt);
result(3,1:2) = r;
ci = nlparci(r,res,'jacobian',J);
result(3,3) = ci(2,2) - r(2);

t1 = 0:max(max(entertime(:,1)),max(exittime(:,1)));

figure(6)
semilogy(entertime(:,1),entertime(:,2),'ob',exittime(:,1),exittime(:,2),'or', ...
    enter2exittime(:,1),enter2exittime(:,2),'og', ...
    t1,efit(result(1,1:2),t1),'-b',t1,efit(result(2,1:2),t1),'-r',...
    t1,efit(result(3,1:2),t1),'-g', ...
    'LineWidth',2.5,'MarkerSize',8);
legend(['Enter: \tau = ' num2str(1/result(1,2),'%1.2f') ' s'], ...
    ['Exit: \tau = ' num2str(1/result(2,2),'%1.2f') ' s'], ...
    ['Enter -> Exit: \tau = ' num2str(1/result(3,2),'%1.2f') ' s']) 
xlabel('time until next event (s)')
ylabel('pdf')
title('distribution of time until next event')
ylim([10^-6 10^0])
set(gca,'FontSize',16);

entercdf = histogram2pts(tenternext,'FaceColor','b','Normalization','cdf','BinMethod','fd');
exitcdf = histogram2pts(texitnext,'FaceColor','r','Normalization','cdf','BinMethod','fd');
enter2exitcdf = histogram2pts(tenterexit,'FaceColor','g','Normalization','cdf','BinMethod','fd');

[r,~,res,~,~,~,J] = lsqcurvefit(efit,[1,1],entercdf(:,1),1-entercdf(:,2),[],[],opt);
resultcdf(1,1:2) = r;
ci = nlparci(r,res,'jacobian',J);
resultcdf(1,3) = ci(2,2) - r(2);

[r,~,res,~,~,~,J] = lsqcurvefit(efit,[1,1],exitcdf(:,1),1-exitcdf(:,2),[],[],opt);
resultcdf(2,1:2) = r;
ci = nlparci(r,res,'jacobian',J);
resultcdf(2,3) = ci(2,2) - r(2);

[r,~,res,~,~,~,J] = lsqcurvefit(efit,[1,1],enter2exitcdf(:,1),1-enter2exitcdf(:,2),[],[],opt);
resultcdf(3,1:2) = r;
ci = nlparci(r,res,'jacobian',J);
resultcdf(3,3) = ci(2,2) - r(2);

figure(7)
semilogy(entercdf(:,1),1-entercdf(:,2),'ob',exitcdf(:,1),1-exitcdf(:,2),'or', ...
    enter2exitcdf(:,1),1-enter2exitcdf(:,2),'og', ...
    t1,efit(resultcdf(1,1:2),t1),'-b',t1,efit(resultcdf(2,1:2),t1),'-r',...
    t1,efit(resultcdf(3,1:2),t1),'-g',...
    'LineWidth',2.5,'MarkerSize',8);
legend(['Enter: \tau = ' num2str(1/resultcdf(1,2),'%1.2f') ' s'], ...
    ['Exit: \tau = ' num2str(1/resultcdf(2,2),'%1.2f') ' s'],...
    ['Enter -> Exit: \tau = ' num2str(1/resultcdf(3,2),'%1.2f') ' s'])  
xlabel('time until next event (s)')
ylabel('1-cdf')
title('time until next event')
ylim([10^-6 10^0])
set(gca,'FontSize',16);

efit2 = @(x,t)x(1).*exp(-x(2).*t)+x(3).*exp(-x(4).*t);

[r,~,res,~,~,~,J] = lsqcurvefit(efit2,[0.5,1,0.5,0.1],entercdf(:,1),1-entercdf(:,2),[],[],opt);
resultcdf2(1,1:4) = r;
ci = nlparci(r,res,'jacobian',J);
resultcdf2(1,5) = ci(2,2) - r(2);
resultcdf2(1,6) = ci(4,2) - r(4);

[r,~,res,~,~,~,J] = lsqcurvefit(efit2,[0.5,1,0.5,0.1],exitcdf(:,1),1-exitcdf(:,2),[],[],opt);
resultcdf2(2,1:4) = r;
ci = nlparci(r,res,'jacobian',J);
resultcdf2(2,5) = ci(2,2) - r(2);
resultcdf2(2,6) = ci(4,2) - r(4);

[r,~,res,~,~,~,J] = lsqcurvefit(efit2,[0.5,1,0.5,0.1],enter2exitcdf(:,1),1-enter2exitcdf(:,2),[],[],opt);
resultcdf2(3,1:4) = r;
ci = nlparci(r,res,'jacobian',J);
resultcdf2(3,5) = ci(2,2) - r(2);
resultcdf2(3,6) = ci(4,2) - r(4);


figure(8)
semilogy(entercdf(:,1),1-entercdf(:,2),'ob',exitcdf(:,1),1-exitcdf(:,2),'or', ...
    enter2exitcdf(:,1),1-enter2exitcdf(:,2),'og', ...
    t1,efit2(resultcdf2(1,1:4),t1),'-b',t1,efit2(resultcdf2(2,1:4),t1),'-r',...
    t1,efit2(resultcdf2(3,1:4),t1),'-g',...
    'LineWidth',2.5,'MarkerSize',8);
legend(['Enter: \tau_1 = ' num2str(1/resultcdf2(1,2),'%1.2f') ' s, \tau_2 = ' num2str(1/resultcdf2(1,4),'%1.2f') ' s, A_1/A_2 = ' num2str(resultcdf2(1,1)/resultcdf2(1,3),'%1.3f')], ...
    ['Exit: \tau_1 = ' num2str(1/resultcdf2(2,2),'%1.2f') ' s, \tau_2 = ' num2str(1/resultcdf2(2,4),'%1.2f') ' s, A_1/A_2 = ' num2str(resultcdf2(2,1)/resultcdf2(2,3),'%1.3f')], ...
    ['Enter -> Exit: \tau_1 = ' num2str(1/resultcdf2(3,2),'%1.2f') ' s, \tau_2 = ' num2str(1/resultcdf2(3,4),'%1.2f') ' s, A_1/A_2 = ' num2str(resultcdf2(3,1)/resultcdf2(3,3),'%1.3f')])  
xlabel('time until next event (s)')
ylabel('1-cdf')
title('time until next event')
ylim([10^-3 10^0])
set(gca,'FontSize',16);
