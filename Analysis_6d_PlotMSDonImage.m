clearvars; 

namestr = '161101_3D3_P11_2b_bottom'; 

singleidx = 11; %which one to plot in it's own figure
singlesize = 300; 
typeclass = 2; %1 for granule; 2 for scrum
sperframe = 0.1; %s per frame usually 0.1s (10 Hz)
figpos = [100 100 600 600];
imgbuffer = 15;
scalebar = 1000; %nm
pixelsize = 71; %nm

lscale = scalebar/pixelsize;

load(['Working/' namestr '_msd_new.mat']); 
load(['Working/' namestr '_sizeinter.mat']); 
I = imread([namestr '_Part1.tif']);

trajclass = ids(2,:) == typeclass; %select the class of trajectories

classtracked = tracked(tracked(:,5) == typeclass,:);
startf1 = classtracked(classtracked(:,3) == 1,4);

nvistraj = length(startf1);

%select the trajectories to make a pretty picture
idxs = zeros(nvistraj,1);
for m=1:nvistraj
    idxs(m) = find((ids(1,:) == startf1(m)) & trajclass);
end
visids = ids(:,idxs);
visids(3,:) = visids(3,:)+1; %add one to offset points vs MSD difference 

coord = zeros(sum(visids(3,:)),3);
piccoord = zeros(nvistraj,2);
startidx = zeros(nvistraj+1,1);  
startidx(2:end) = cumsum(visids(3,:));
startidx = startidx+1;
imgs = zeros([2*imgbuffer+1 2*imgbuffer+1 nvistraj],'like',I);

%collect all the trajectories and subtract each's center of mass to 
%move to the same cordinate system 
for m=1:nvistraj
    %select the correct coordinates from tracked
    subtracked = tracked(tracked(:,4) == visids(1,m),[1,2,4,5]);
    subtracked = subtracked(subtracked(:,4) == visids(2,m),1:3);
    
    %avgx = sum(subtracked(:,1)) / visids(3,m);
    %avgy = sum(subtracked(:,2)) / visids(3,m);
    
    piccoord(m,:) = round([subtracked(1,1),subtracked(1,2)]); 
    imgs(:,:,m) = ...
        I(piccoord(m,2)-imgbuffer:piccoord(m,2)+imgbuffer,piccoord(m,1)-imgbuffer:piccoord(m,1)+imgbuffer);
    
    coord(startidx(m):startidx(m+1)-1,1) = subtracked(:,3);
    coord(startidx(m):startidx(m+1)-1,2) = subtracked(:,1)-subtracked(1,1);
    coord(startidx(m):startidx(m+1)-1,3) = subtracked(:,2)-subtracked(1,2);
end


col = 0:sperframe:max(visids(3,:))*sperframe;

% big = max(max(coord(:,2:3)));
% small = min(min(coord(:,2:3)));
% limit = max(abs([big small]));
ncolors = 64;
custommap = [gray(ncolors);parula(ncolors)];

if ishghandle(1)
    close(1)
end
figure(1)
set(gcf,'Position',figpos,'Color',[1 1 1]);
gridn = ceil(sqrt(nvistraj));
mfct = 1/gridn;
for m=1:nvistraj
    curr = coord(coord(:,1) == visids(1,m),:);
    big = max(max(curr(:,2:3)));
    small = min(min(curr(:,2:3)));
    limit = max(abs([big small]))+2;
    if limit < imgbuffer
        limit = imgbuffer;
    end
    z = zeros(size(curr(:,2)))';
    xs = mfct.*mod(m-1,gridn);
    ys = mfct.*floor((m-1)/gridn);
    h(m) = subplot('Position',[xs+0.001 1-mfct-ys+0.001 mfct-0.002 mfct-0.002]);
    colormap(custommap);
    ph(1) = imagesc([-imgbuffer imgbuffer],[-imgbuffer imgbuffer],imgs(:,:,m));
    hold on
    ph(2) = surface([curr(:,2)';curr(:,2)'], ...
        [curr(:,3)';curr(:,3)'], ...
        [z;z],[col(1:length(curr));col(1:length(curr))],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',3);
    hold off
    set(gca,'Color',[0 0 0],'xtick',[],'ytick',[],'YDir','normal');
    oC1 = double(get(ph(1),'CData'));
    cmin1 = min(min(oC1));
    cmax1 = max(max(oC1));
    C1 = min(ncolors,round((ncolors-1)*(oC1-cmin1)/(cmax1-cmin1))+1);
    set(ph(1),'CData',C1);
    oC2 = get(ph(2),'CData');
    cmin2 = min(min(oC2));
    cmax2 = max(max(oC2));
    C2 = min(ncolors,round((ncolors-1)*(oC2-cmin2)/(cmax2-cmin2))+1);
    C2 = C2+ncolors;
    set(ph(2),'CData',C2);
    xlim([-limit limit])
    ylim([-limit limit])
    caxis([min(min(C1)) max(max(C2))])
    axis square
    line([limit-1-lscale limit-1],[-limit+3 -limit+3],'Color','w','LineWidth',6)
 


%     text(-limit+2,limit-2,[num2str(m) '. ' num2str(visids(1,m))],'Color','w','FontSize',14);
%     if(m == 1)
%         ostr = strrep(namestr,'_',' ');
%         text(-limit,-limit+2,[ostr ': ' num2str(vistrajl)],'Color','w','FontSize',12);
%     end
end
tightfig;



if ishghandle(2)
    close(2)
end
figure(2)
set(gcf,'Position',[100 100 singlesize singlesize],'Color',[1 1 1]);

curr = coord(coord(:,1) == visids(1,singleidx),:);
big = max(max(curr(:,2:3)));
small = min(min(curr(:,2:3)));
limit = max(abs([big small]))+2;
if limit < imgbuffer
    limit = imgbuffer;
end
z = zeros(size(curr(:,2)))';
colormap(custommap);
ph(1) = imagesc([-imgbuffer imgbuffer],[-imgbuffer imgbuffer],imgs(:,:,m));
hold on
ph(2) = surface([curr(:,2)';curr(:,2)'], ...
    [curr(:,3)';curr(:,3)'], ...
    [z;z],[col(1:length(curr));col(1:length(curr))],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3);
hold off
set(gca,'Color',[0 0 0],'xtick',[],'ytick',[],'YDir','normal');
oC1 = double(get(ph(1),'CData'));
cmin1 = min(min(oC1));
cmax1 = max(max(oC1));
C1 = min(ncolors-1,round((ncolors-1)*(oC1-cmin1)/(cmax1-cmin1))+1);
set(ph(1),'CData',C1);
oC2 = get(ph(2),'CData');
cmin2 = min(min(oC2));
cmax2 = max(max(oC2));
C2 = min(ncolors,round((ncolors-1)*(oC2-cmin2)/(cmax2-cmin2))+1);
C2 = C2+ncolors;
set(ph(2),'CData',C2);
xlim([-limit limit])
ylim([-limit limit])
caxis([min(min(C1)) max(max(C2))])
axis square
line([limit-1-lscale limit-1],[-limit+3 -limit+3],'Color','w','LineWidth',6)

tightfig;

