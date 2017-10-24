clearvars; 

namestr = 'Simulated1'; 
fitlength = 10; %number of points to fit for MSD, also min MSD to be included
interactioncutoff = 22; %pixel difference to be considered interacting 
pixelsize = 71; %nm from microscope
sperframe = 0.1; %s per frame usually 0.1s (10 Hz)

load(['Working/' namestr '_sizeinter.mat']); 
load(['Working/' namestr '_msd_new.mat']); 

%filter the trajectories for the required size 
idxs = find(ids(3,:) >= fitlength); 
ntraj = length(idxs);
nsmall = sum(ids(2,idxs) == 1);
nlarge = sum(ids(2,idxs) == 2);
disp(['Found ' num2str(ntraj) ' trajectories of at least length ' num2str(fitlength) '.']);
disp(['Of which small: ' num2str(nsmall) ' large: ' num2str(nlarge)]);
disp(['Excluded ' num2str(size(ids,2)-1-ntraj) ' trajectories.']);

%create variables to go into a table to export to R 
ID = ids(1,idxs)'; %done
sizetype = char(zeros(ntraj,1)); %edit in loop
D = zeros(ntraj,1);
alpha = zeros(ntraj,1);
tlength = ids(3,idxs)'; %done
gsize = zeros(ntraj,1); %calc in loop
sizeerrors = zeros(ntraj,1); %calc in loop 
smalldist = zeros(ntraj,1); %calc in loop
interacting = zeros(ntraj,1); %calc in loop 
t = sperframe.*(1:fitlength);
logt = log10(t);
opt = optimset('Display','none','TolFun',10^-8,'TolX',10^-8);
lfit = @(x,t)x(1).*t+x(2); 
wh = waitbar(0,'Starting...'); 
tic;
ud = 100;
for m=1:ntraj
    %calc the other variables for the table first 
    if ids(2,idxs(m)) == 1
        sizetype(m) = 'S';
    else
        sizetype(m) = 'L';
    end
    pts = tracked(and(tracked(:,4) == ids(1,idxs(m)),tracked(:,5) == ids(2,idxs(m))),6:8);
    npts = size(pts,1);
    smalldist(m) = min(pts(:,1));
    pts(:,4) = double(pts(:,1) <= interactioncutoff);
    spts = squeeze(sum(pts,1));
    gsize(m) = spts(2)/npts;
    sizeerrors(m) = spts(3);
    interacting(m) = spts(4);
    
    %now fit the msd
    cmsd = msd(1:fitlength,idxs(m))';
    logmsd = log10(cmsd);
    result = lsqcurvefit(lfit,[1,1],logt,logmsd,[],[],opt);
    alpha(m) = result(1);
    D(m) = result(2); %in units of log10(pixel2/s)
    if ud == m
        ctime = toc;
        tperitr = ctime/m;
        tleft = round(tperitr*(ntraj-m));
        hr = floor(tleft/3600);
        mins = floor((tleft-hr*3600)/60);
        sec = round(tleft-hr*3600-mins*60);
        waitbar(m/ntraj,wh,{['Fitting Trajectory: ' num2str(m) ' of ' num2str(ntraj)]; ...
            ['Time Remaining: ' num2str(hr,'%02i') ':' num2str(mins,'%02i') ':' num2str(sec,'%02i')]});
        ud = ud + 100;
    end
end

%convert size from pixels to distance & std to FWHM
cFWHM = 2*sqrt(2*log(2));
smalldist = pixelsize.*smalldist;
gsize = cFWHM*pixelsize.*gsize;

%convert D to real units
D2 = 10.^D; %units of pixel^2/s
D2 = D2.*pixelsize^2*(1/1000)^2/4; %now um^2/s + divide by 4 for 2-D system 

%create table 
T = table(ID,sizetype,D2,alpha,tlength,gsize,sizeerrors,smalldist,interacting, ...
    'VariableNames',{'ID','SizeClass','D','alpha','Length','Size','NumSizeErrors','SmallestDistance','NumFramesInteracting'});
Tfull = table(ids(1,2:end)',ids(2,2:end)',ids(3,2:end)', ...
    'VariableNames',{'ID','SizeClass','Length'});


writetable(T,['Output/' namestr '_fit' num2str(fitlength) '_inter' num2str(interactioncutoff) '__new_R.csv']);
writetable(Tfull,['Output/' namestr '_fulllength.csv']);
save(['Working/' namestr '_fit' num2str(fitlength) '_inter' num2str(interactioncutoff) '_fulldatanew.mat'], ...
    'T','idxs','tracked','msd','ids','fitlength','namestr','interactioncutoff','D');
%use the above save for future GUI to explore results. 

close(wh);