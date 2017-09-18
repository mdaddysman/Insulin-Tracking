clearvars; 

namestr = '161101_3D3_P11_2b_bottom'; 

vistrajid = [25969, 59137, 38034, 60977]; %161101_3D3_P11_2b_bottom
%vistrajid = [29153, 18693, 19670, 9470]; %161101_3D3_P11_2c_top
%vistrajid = [27888, 29872, 25208, 38562]; %161101_3D3_P11_2d
typeclass = 1;
sperframe = 0.1; %s per frame usually 0.1s (10 Hz)

nvistraj = length(vistrajid);

load(['Working/' namestr '_msd_new.mat']); 
load(['Working/' namestr '_sizeinter.mat']); 

trajclass = ids(2,:) == typeclass; %select the class of trajectories 

%select the trajectories to make a pretty picture
idxs = zeros(nvistraj,1);
for m=1:nvistraj
    idxs(m) = find((ids(1,:) == vistrajid(m)) & trajclass);
end
visids = ids(:,idxs);
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


col = 0:sperframe:max(visids(3,:))*sperframe;

big = max(max(coord(:,2:3)));
small = min(min(coord(:,2:3)));
limit = max(abs([big small]));

if ishghandle(3)
    close(3)
end
figure(3)
set(gcf,'Position',[100 100 400 400],'Color',[1 1 1]);
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
        'linew',3);
    set(gca,'Color',[0 0 0],'xtick',[],'ytick',[],'CLim',[0 max(visids(3,:))*sperframe]);
    xlim([-limit limit])
    ylim([-limit limit])
    axis square
    if(m == nvistraj)
        line([limit-1-14.0845 limit-1],[-limit+3 -limit+3],'Color','w','LineWidth',6)
    end
%     text(-limit+2,limit-2,[num2str(m) '. ' num2str(visids(1,m))],'Color','w','FontSize',14);
%     if(m == 1)
%         ostr = strrep(namestr,'_',' ');
%         text(-limit,-limit+2,[ostr ': ' num2str(vistrajl)],'Color','w','FontSize',12);
%     end
end
tightfig;