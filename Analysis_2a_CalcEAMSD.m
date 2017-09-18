clearvars;

namestr = '161101_3D3_P11_2d';

load(['Working/' namestr '_sizeinter.mat']); 

sperframe = 0.1; 

%eamsd gives the following columns
%time small s=1,2,3,4 large s=1,2,3,4 for 9 columns 
%will calculate for max time - 5 time steps from the end 

%set up inital variables 
nsmall = max(tracked(tracked(:,5) == 1,4));
nlarge = max(tracked(tracked(:,5) == 2,4));
ntraj = nsmall + nlarge;
nframes = max(tracked(:,3));

eamsd = zeros([nframes-5 9]);
eamsd(:,1) = sperframe.*(0:nframes-6); 
eatrmsd = zeros([nframes-1 9]);
eatrmsd(:,1) = sperframe.*(0:nframes-2); 
msdgrid = zeros([nframes nsmall 4]); 
msdtrgrid = zeros([nframes nsmall 4]); 

wh = waitbar(0,'Starting...'); 
tic;
subpts = tracked(tracked(:,5) == 1,:);
nshift = 0; 
sorl = 1; %small or large
ud = 100;
for n=1:ntraj
    tpts = subpts(subpts(:,4) == n-nshift,1:3);
    tpts = sortrows(tpts,3);
    for m=1:4
        clength = size(tpts,1);
        if clength <= m
            continue;
        end
       temp = (tpts(1+m:end,1) - tpts(1:end-m,1)).^2 + (tpts(1+m:end,2) - tpts(1:end-m,2)).^2;
       
       msdgrid(tpts(1,3):tpts(1,3)+length(temp)-1,n-nshift,m) = temp;  
       msdtrgrid(1:length(temp),n-nshift,m) = temp;
    end
    
    if ud == n
        ctime = toc;
        tperitr = ctime/n;
        tleft = round(tperitr*(ntraj-n));
        hr = floor(tleft/3600);
        mins = floor((tleft-hr*3600)/60);
        sec = round(tleft-hr*3600-mins*60);
        waitbar(n/ntraj,wh,{['MSD Trajectory: ' num2str(n) ' of ' num2str(ntraj)]; ...
            ['Time Remaining: ' num2str(hr,'%02i') ':' num2str(mins,'%02i') ':' num2str(sec,'%02i')]});
        ud = ud + 100;
    end
    if n == nsmall %move onto large
        lmat = msdgrid > 0;
        ltrmat = msdtrgrid > 0;
        for z=1:4
            summsd = squeeze(sum(msdgrid(:,:,z),2));
            sumnorm = squeeze(sum(lmat(:,:,z),2));
            eamsd(:,z+1) = summsd(1:end-5)./sumnorm(1:end-5);
            sumtrmsd = squeeze(sum(msdtrgrid(:,:,z),2));
            sumtrnorm = squeeze(sum(ltrmat(:,:,z),2));
            eatrmsd(:,z+1) = sumtrmsd(1:end-1)./sumtrnorm(1:end-1);
        end
        subpts = tracked(tracked(:,5) == 2,:);
        nshift = nsmall;
        msdgrid = zeros([nframes nlarge 4]); 
        msdtrgrid = zeros([nframes nlarge 4]);
        sorl = 2;
    end
end

lmat = msdgrid > 0;
ltrmat = msdtrgrid > 0;
for z=1:4
    summsd = squeeze(sum(msdgrid(:,:,z),2));
    sumnorm = squeeze(sum(lmat(:,:,z),2));
    eamsd(:,z+5) = summsd(1:end-5)./sumnorm(1:end-5);
    sumtrmsd = squeeze(sum(msdtrgrid(:,:,z),2));
    sumtrnorm = squeeze(sum(ltrmat(:,:,z),2));
    eatrmsd(:,z+5) = sumtrmsd(1:end-1)./sumtrnorm(1:end-1);
    
end

save(['Working/' namestr '_eamsd.mat'],'eamsd','eatrmsd'); 
close(wh);