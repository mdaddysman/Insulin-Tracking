namestr = '170427_3B11M_P13_Plate2a_Top';
nparts = 3; 

for m=1:nparts
    d = csvread(['Input Tracking/' namestr '_Part' num2str(m) '_s_Results.csv']);
    d = d(:,3:12);
    d(:,1) = d(:,1)+1;
    if(m == 1)
        alldata = d;
    else
        d(:,1) = d(:,1) + max(alldata(:,1));
        alldata = [alldata; d];
    end
end

for m=1:nparts
    d = csvread(['Input Tracking/' namestr '_Part' num2str(m) '_l_Results.csv']);
    d = d(:,3:12);
    d(:,1) = d(:,1)+1;
    if(m == 1)
        largedata = d;
    else
        d(:,1) = d(:,1) + max(largedata(:,1));
        largedata = [largedata; d];
    end
end

save(['Working/' namestr '_combine.mat'],'alldata','largedata'); 