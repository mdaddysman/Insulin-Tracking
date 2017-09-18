namestr = '170427_3B11M_P13_Plate2a_Top';
linkrange = 9; 

load(['Working/' namestr '_SPIFF.mat']); 

%prepare the data for tracking 
smallpos = sortrows(smalldata(:,[2:3 1]),3);
largepos = sortrows(largedata(:,[2:3 1]),3);

trackedsmall = track(smallpos,linkrange);
trackedlarge = track(largepos,linkrange);

save(['Working/' namestr '_traj.mat'],'trackedlarge','trackedsmall','linkrange');

%Move onto Analysis