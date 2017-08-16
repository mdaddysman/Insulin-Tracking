clearvars;

namestr = '170427_3B11M_P13_Plate2a_Top';

load([namestr '_sizeinter.mat']); 

%ids give the trajectory number, large or small, length 
%msd matches the index of ids col 1 is the time displacement in frame #
%therefore col 1 in ids is NaN 

%set up inital variables 
nsmall = max(tracked(tracked(:,5) == 1,4));
nlarge = max(tracked(tracked(:,5) == 2,4));
ntraj = nsmall + nlarge;
nframes = max(tracked(:,3));
ids = zeros([3 ntraj+1]); 
msd = zeros([nframes-1 ntraj+1]);
msd(:,1) = 1:nframes-1; 
ids(1:3,1) = [NaN NaN NaN];

wh = waitbar(0,'Starting...'); 
tic;
subpts = tracked(tracked(:,5) == 1,:);
nshift = 0; 
sorl = 1; %small or large
ud = 100;
for n=1:ntraj
    tpts = subpts(subpts(:,4) == n-nshift,1:3);
    ids(1,n+1) = n-nshift;
    ids(2,n+1) = sorl;
    ids(3,n+1) = size(tpts,1)-1;
    tpts = sortrows(tpts,3);
    %sanity check
    delframe = tpts(2:end,3) - tpts(1:end-1,3);
    if sum(delframe) ~= ids(3,n+1)
        disp(['Trajectory ' num2str(n) ' is invalid.']);
    end
    for m=1:ids(3,n+1)
        idx1 = 1:m+1:ids(3,n+1);
        idx2 = m+1:m+1:ids(3,n+1);
        if length(idx1) > length(idx2)
            idx1 = idx1(1:end-1);
        end
        if length(idx2) > length(idx1)
            idx2 = idx2(1:end-1);
        end
        if (isempty(idx1) && isempty(idx2))
            idx1 = 1;
            idx2 = ids(3,n+1);
        end
        temp = (tpts(idx2,1) - tpts(idx1,1)).^2 + (tpts(idx2,2) - tpts(idx1,2)).^2;
        msd(m,n+1) = sum(temp)/length(temp);
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
        subpts = tracked(tracked(:,5) == 2,:);
        nshift = nsmall;
        sorl = 2;
    end
end

save([namestr '_msd_new.mat'],'msd','ids'); 
close(wh);