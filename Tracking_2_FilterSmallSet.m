namestr = '170427_3B11M_P13_Plate2a_Top';
radius = 10;  

load([namestr '_combine.mat']);

nframes = max(alldata(:,1));

wh = waitbar(0,['Processing frame 0 of ' num2str(nframes)]);

for m=1:nframes %for each frame 
    waitbar(m/nframes,wh,['Processing frame ' num2str(m) ' of ' num2str(nframes)]);
    bigcenters = largedata(m == largedata(:,1),2:3); %pull out the centers of the circles in each frame
    smallidx = find(alldata(:,1) == m);
    fspts = alldata(smallidx,2:3);
    exclude = false(length(smallidx),1); %keep track if the point is inside one of the circles 
    for n=1:size(bigcenters,1) %now go through all of the big centers in the frame 
        dist = sqrt((fspts(:,1) - bigcenters(n,1)).^2   +   (fspts(:,2) - bigcenters(n,2)).^2);
        insr = dist<radius;
        exclude = exclude | insr;
    end
    
    alldata(smallidx,4) = exclude;
        
end

smalldata = alldata(alldata(:,4) == 0,:);

save([namestr '_filtered.mat'],'smalldata','largedata','alldata');
close(wh);
    
    
    