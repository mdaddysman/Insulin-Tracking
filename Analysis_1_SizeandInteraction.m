clearvars;
global gridX gridY fitpicture
namestr = '170427_3B11M_P13_Plate2a_Top';
nparts = 3;
pixelsize = 71; %nm - set by the microscope
boxsize = 30; %in pixels (should be an even number)
ctol = 2; %the +/- tolerence for the center to be considered a valid result

runsizeanalysis = false;

%this script adds 5 columns & combines the tracked variables
%col 5 is a small (1) or large (2) particle trajectory
%col 6 is the distance to the closest particle in pixel units
%col 7 is the size of the particle from a gaussian fit
%col 8 is a boolean if the size is "valid" (near the center of the image
%box)
%col 9 is the ID for looking up the fit results in the viewer (1a)

load(['Working/' namestr '_traj.mat']);

%expand the size of the matrix
trackedsmall(:,5:9) = zeros([size(trackedsmall,1) 5]);
trackedlarge(:,5:9) = zeros([size(trackedlarge,1) 5]);

trackedsmall(:,5) = ones([size(trackedsmall,1) 1]);
trackedlarge(:,5) = 2.*ones([size(trackedlarge,1) 1]);

tracked = [trackedsmall; trackedlarge];

%sort by frame - restore before saving
tracked = sortrows(tracked,3);

%find nearest neighbors
wh = waitbar(0,'Distance Analysis Starting...');
nframes = max(tracked(:,3));
tic;

for m=1:nframes
    fidxs = find(tracked(:,3) == m);
    tempdist = zeros(length(fidxs),1);
    fpts = tracked(fidxs,1:2);
    for n=1:length(fidxs)
        dist = sqrt((fpts(n,1)-fpts(:,1)).^2 + (fpts(n,2)-fpts(:,2)).^2);
        dist(n) = [];
        tempdist(n) = min(dist);
    end
    tracked(fidxs,6) = tempdist;
    ctime = toc;
    tperitr = ctime/m;
    tleft = round(tperitr*(nframes-m));
    hr = floor(tleft/3600);
    mins = floor((tleft-hr*3600)/60);
    sec = round(tleft-hr*3600-mins*60);
    waitbar(m/nframes,wh,{['Distance Analysis Frame: ' num2str(m) ' of ' num2str(nframes)]; ...
        ['Time Remaining: ' num2str(hr,'%02i') ':' num2str(mins,'%02i') ':' num2str(sec,'%02i')]});
end

if runsizeanalysis
    
    waitbar(0,wh,'Size Analysis: Creating Point Image Stack');
    
    Img = LoadTiffStack([namestr '_Part1.tif']);
    nimg = size(Img,3);
    
    npoints = size(tracked,1);
    ptImg = zeros([boxsize+1 boxsize+1 npoints],'like',Img);
    corners = zeros(npoints,3);
    %pull out each stack of points
    for m=1:npoints
        corners(m,1) = round(tracked(m,1)-boxsize/2);
        if corners(m,1)+boxsize > size(Img,2)
            corners(m,1) = size(Img,2)-boxsize;
        end
        if corners(m,1) < 1
            corners(m,1) = 1;
        end
        corners(m,2) = round(tracked(m,2)-boxsize/2);
        if corners(m,2)+boxsize > size(Img,1)
            corners(m,2) = size(Img,1)-boxsize;
        end
        if corners(m,2) < 1
            corners(m,2) = 1;
        end
        corners(m,3) = tracked(m,3);
        %     if size(temp,1) < boxsize+1 || size(temp,2) < boxsize+1
        %         temp2 = zeros([boxsize+1 boxsize+1],'like',Img);
        %         temp2(1:size(temp,1),1:size(temp,2)) = temp;
        %         temp = temp2;
        %     end
        
    end
    
    for n=1:nparts
        if n ~= 1
            Img = LoadTiffStack([namestr '_Part' num2str(n) '.tif']);
            nimg(n) = size(Img,3);
            pframe = sum(nimg(1:n-1));
            mstart = mend+1;
        else
            pframe = 0;
            mstart = 1;
        end
        mend = find(corners(:,3) == sum(nimg),1,'last');
        
        for m=mstart:mend
            ptImg(:,:,m) = Img(corners(m,2):corners(m,2)+boxsize,corners(m,1):corners(m,1)+boxsize,corners(m,3)-pframe);
        end
        clear Img
    end
    
    %fit 2D gaussian
    [gridX,gridY] = meshgrid(1:boxsize+1);
    opt = optimset('Display','none','TolFun',10^-8,'TolX',10^-8);  % 'none' or 'iter'
    t = 0:pixelsize:boxsize*pixelsize;
    tfit = 0:pixelsize/4:boxsize*pixelsize;
    pfit = 1:0.25:boxsize+1;
    results = zeros([12 npoints],'double');
    tolmax = boxsize/2 + ctol;
    tolmin = boxsize/2 - ctol;
    ud = 100;
    tic;
    for m=1:npoints
        fitpicture = double(ptImg(:,:,m));
        if tracked(m,5) == 2
            wguess = 8;
        else
            wguess = 3;
        end
        [result,resnorm] = lsqnonlin('gauss2D',[max(max(fitpicture)),boxsize/2,boxsize/2,wguess,wguess,100], ...
            [0,0,0,0,0,0],[70000,boxsize,boxsize,2*boxsize,2*boxsize,70000],opt);
        results(1,m) = m;
        results(2:7,m) = result;
        results(8,m) = resnorm;
        results(9:11,m) = tracked(m,3:5); %frame, traj #, large or small
        
        tracked(m,7) = sum(result(4:5))/2;
        if result(2) <= tolmax && result(2) >= tolmin && result(3) <= tolmax && result(3) >= tolmin
            tracked(m,8) = 1;
            results(12,m) = 1;
        end
        tracked(m,9) = m;
        
        if(m == ud)
            ctime = toc;
            tperitr = ctime/m;
            tleft = round(tperitr*(npoints-m));
            hr = floor(tleft/3600);
            mins = floor((tleft-hr*3600)/60);
            sec = round(tleft-hr*3600-mins*60);
            waitbar(m/npoints,wh,{['Size Fitting Point: ' num2str(m) ' of ' num2str(npoints)]; ...
                ['Time Remaining: ' num2str(hr,'%02i') ':' num2str(mins,'%02i') ':' num2str(sec,'%02i')]});
            ud = ud+100;
        end
    end
    results = results';
end
close(wh);
%restore order - sort by traj, frame, large small
tracked = sortrows(tracked,[4 3 5]);
save(['Working/' namestr '_sizeinter.mat'],'linkrange','tracked','runsizeanalysis');
if runsizeanalysis
    save(['Working/' namestr '_fits.mat'],'results','t','tfit','boxsize','pixelsize','ptImg');
end