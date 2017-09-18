namestr = 'cMovie1';

rs = 4;
rl = 10; 
nframes = 100; 
width = 860; height = 800; 

load(['Working/' namestr '_filtered.mat']); 

image = zeros([height width 3 nframes],'uint8'); 

excludedata = alldata(alldata(:,4) == 1,1:3);
[yy,xx] = meshgrid(1:width,1:height);

wh = waitbar(0,['Processing frame 0 of ' num2str(nframes)]);
for m=1:nframes
    waitbar(m/nframes,wh,['Processing frame ' num2str(m) ' of ' num2str(nframes)]);
    smallc = round(smalldata(smalldata(:,1) == m,2:3));
    largec = round(largedata(largedata(:,1) == m,2:3)); 
    excludec = round(excludedata(excludedata(:,1) == m,2:3));
    
    %now go through each center and make the circle 
    Cf = false(height,width);
    for n=1:size(smallc,1)  
        C = sqrt((yy-smallc(n,2)).^2 + (xx-smallc(n,1)).^2)==rs; 
        Cf = Cf | C;
    end
    image(:,:,1,m) = uint8(255.*Cf);  %plot small circles in red
    
    Cf = false(height,width);
    for n=1:size(largec,1)
        C = sqrt((yy-largec(n,2)).^2 + (xx-largec(n,1)).^2)==rl;
        Cf = Cf | C;
    end
    image(:,:,2,m) = uint8(255.*Cf);  %plot large circles in green
    
    Cf = false(height,width);
    for n=1:size(excludec,1)
        C = sqrt((yy-excludec(n,2)).^2 + (xx-excludec(n,1)).^2)==rs;
        Cf = Cf | C;
    end
    image(:,:,3,m) = uint8(255.*Cf);  %plot small excluded circles in blue
    
end

DisplayStack(image)
close(wh)