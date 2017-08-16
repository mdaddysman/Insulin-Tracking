outputname = '161101_3D3_P11_middle_granule.csv';

data = trackedsmall; 
secperframe = 0.1; %seconds per frame of movie 

ndata = zeros(size(data,1),5);

ndata(:,1) = data(:,4); %traj ID
ndata(:,2) = secperframe.*(data(:,3)-1); %absolute time
ndata(:,4:5) = data(:,1:2); %X,Y positions 

for m=1:max(ndata(:,1)) %calc relative time 
    idxs = find(ndata(:,1) == m);
    
    ndata(idxs,3) = ndata(idxs,2) - ndata(idxs(1),2);
end


csvwrite(outputname,ndata);