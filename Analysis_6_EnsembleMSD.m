clearvars; 

namestr = '161101_3D3_P11_2b_bottom'; 
fitlength = 30; %number of points to fit for MSD, 
minlength = 10; %min MSD to be included
nvistraj = 16;
vistrajl = 200;
pixelsize = 71; %nm from microscope
sperframe = 0.1; %s per frame usually 0.1s (10 Hz)
xlimit = [0 10];

load(['Working/' namestr '_msd_new.mat']); 
load(['Working/' namestr '_sizeinter.mat']); 

trajclass = ids(2,:) == 1; %select the class of trajectories 
trajlength = ids(3,:) >= minlength; %select the long enough trajectories 
trajboth = trajclass & trajlength; %combine

idxboth = find(trajboth);
idxclass = find(trajclass); 


msdboth = msd(:,idxboth);
msdclass = msd(:,idxclass);

cpixel = pixelsize^2*(1/1000)^2; %pixel conversion to um^2/s

avgmsdboth = sum(msdboth,2) ./ sum(msdboth>0,2) .* cpixel;
avgmsdclass = sum(msdclass,2) ./ sum(msdclass>0,2) .* cpixel;


t = 1:size(msd,1)/2;
t = sperframe.*t; %convert from frames to s 

%fit the MSD
logt = log10(t(1:fitlength));
opt = optimset('Display','none','TolFun',10^-8,'TolX',10^-8);
lfit = @(x,t)x(1).*t+x(2);
logmsd = log10(avgmsdboth(1:fitlength));
[resultboth,resnorm,res,~,~,~,J] = lsqcurvefit(lfit,[1,1],logt,logmsd',[],[],opt);
ci = nlparci(resultboth,res,'jacobian',J);
plusminusboth(1,:) = ci(1,:) - resultboth(1);
plusminusboth(2,:) = ci(2,:) - resultboth(2);


logmsd = log10(avgmsdclass(1:fitlength));
resultclass = lsqcurvefit(lfit,[1,1],logt,logmsd',[],[],opt);

fitboth = 10^resultboth(2).*t.^resultboth(1);
fitclass = 10^resultclass(2).*t.^resultclass(1);

figure(1)
subplot(2,1,1)
plot(t,avgmsdboth(1:length(t)),'o',t,fitboth,'-')
title(['Granules and longer than ' num2str(fitlength) ' frames'])
xlim(xlimit)
xlabel('delta (s)')
ylabel('MSD (um^2)')

subplot(2,1,2)
plot(t,avgmsdclass(1:length(t)),'o',t,fitclass,'-')
title('Granules of all trajectory lengths') 
xlim(xlimit)
xlabel('delta (s)')
ylabel('MSD (um^2)')

save(['Working/' namestr '_class' num2str(ids(2,idxclass(1))) '_ensemblefit.mat'], ...
    't','avgmsdboth','avgmsdclass','fitboth','fitclass','resultboth','resultclass','fitlength','minlength','plusminusboth');

%select the trajectories to make a pretty picture 
vtrajlength = ids(3,:) >= vistrajl;
idxs = find(vtrajlength & trajclass);
shuffle = randperm(length(idxs),nvistraj); 
nidxs = idxs(shuffle); 
visids = ids(:,nidxs);
visids(3,:) = visids(3,:)+1; %add one to offset points vs MSD difference 

coord = zeros(sum(visids(3,:)),3);
startidx = zeros(nvistraj+1,1);  
startidx(2:end) = cumsum(visids(3,:));
startidx = startidx+1;

%collect all the trajectories and subtract each's center of mass to 
%move to the same cordinate system 
for m=1:nvistraj
    %select the correct coordinates from tracked
    subtracked = tracked(tracked(:,4) == visids(1,m),[1,2,4,5]);
    subtracked = subtracked(subtracked(:,4) == visids(2,m),1:3);
    
    avgx = sum(subtracked(:,1)) / visids(3,m);
    avgy = sum(subtracked(:,2)) / visids(3,m);
    
    coord(startidx(m):startidx(m+1)-1,1) = subtracked(:,3);
    coord(startidx(m):startidx(m+1)-1,2) = subtracked(:,1)-avgx;
    coord(startidx(m):startidx(m+1)-1,3) = subtracked(:,2)-avgy;
end

% coord2 = coord;
% minx = min(coord(:,2));
% if minx < 1
%     coord2(:,2) = coord(:,2) + abs(minx) + 1;
% end
% miny = min(coord(:,3));
% if miny < 1
%     coord2(:,3) = coord(:,3) + abs(miny) + 1;
% end
% 
% rcoord = coord2;
% rcoord(:,2:3) = round(rcoord(:,2:3).*10);
% 
% largestd = ceil(max(max(rcoord(:,2:3))));
% 
% trajimgstack = zeros([largestd largestd 3 nvistraj],'uint8');
% 
% %add color scale here later 
% cmap = parula(max(visids(3,:)));
% cmap = 255.*cmap; 
% 
% for m=1:nvistraj
%     curr = rcoord(rcoord(:,1) == visids(1,m),:);
%     for n=1:length(curr)
%         trajimgstack(curr(n,3),curr(n,2),1:3,m) = cmap(n,1:3);
%     end
% end
% 
% for m=1:3
%     tileimg(:,:,m) = [trajimgstack(:,:,m,1),trajimgstack(:,:,m,2),trajimgstack(:,:,m,3),trajimgstack(:,:,m,4); ...
%         trajimgstack(:,:,m,5),trajimgstack(:,:,m,6),trajimgstack(:,:,m,7),trajimgstack(:,:,m,8); ...
%         trajimgstack(:,:,m,9),trajimgstack(:,:,m,10),trajimgstack(:,:,m,11),trajimgstack(:,:,m,12); ...
%         trajimgstack(:,:,m,13),trajimgstack(:,:,m,14),trajimgstack(:,:,m,15),trajimgstack(:,:,m,16)];
% end
% 
% figure(2)
% imshow(tileimg)

col = 0:sperframe:max(visids(3,:))*sperframe;

big = max(max(coord(:,2:3)));
small = min(min(coord(:,2:3)));
limit = max(abs([big small]));

if ishghandle(3)
    close(3)
end
figure(3)
set(gcf,'Position',[100 100 1200 1200],'Color',[1 1 1]);
gridn = ceil(sqrt(nvistraj));
mfct = 1/gridn;
for m=1:nvistraj
    curr = coord(coord(:,1) == visids(1,m),:);
    z = zeros(size(curr(:,2)))';
    xs = mfct.*mod(m-1,gridn);
    ys = mfct.*floor((m-1)/gridn);
    h = subplot('Position',[xs+0.001 1-mfct-ys+0.001 mfct-0.002 mfct-0.002]);
    surface([curr(:,2)';curr(:,2)'], ...
        [curr(:,3)';curr(:,3)'], ...
        [z;z],[col(1:length(curr));col(1:length(curr))],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    set(gca,'Color',[0 0 0],'xtick',[],'ytick',[],'CLim',[0 max(visids(3,:))*sperframe]);
    xlim([-limit limit])
    ylim([-limit limit])
    axis square
    text(-limit+2,limit-2,[num2str(m) '. ' num2str(visids(1,m))],'Color','w','FontSize',14);
    if(m == 1)
        ostr = strrep(namestr,'_',' ');
        text(-limit,-limit+2,[ostr ': ' num2str(vistrajl)],'Color','w','FontSize',12);
    end
    if(m == nvistraj)
        line([limit-1-14.0845 limit-1],[-limit+3 -limit+3],'Color','w','LineWidth',6)
    end
end
tightfig;