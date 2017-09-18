clearvars;

namestr = 'cMovie1';

validitylength = 10; %how many frames must the scrum or granule last to be included in the analysis
radius = 17; %how many pixels must the granule be in to be considered interacting 

load(['Working/' namestr '_sizeinter.mat']); 

scrums = tracked(tracked(:,5) == 2,1:4);
granules = tracked(tracked(:,5) == 1,1:4);

clear tracked

vgranules = zeros(max(granules(:,4)),2); 
vgranules(:,1) = 1:max(granules(:,4));

wh = waitbar(0,'Analyzing Granules');
for m=1:max(granules(:,4))
    temp = granules(granules(:,4) == m,:);
    if size(temp,1) >= validitylength
        vgranules(m,2) = 1;
    end
end

vgranules = vgranules(vgranules(:,2) == 1,1);
nvgranules = length(vgranules);

waitbar(0,wh,['Analyzing Granules Found: ' num2str(nvgranules)])

potinter = zeros(nvgranules*2,9); %matrix of potential interactions
%1-Granule#, 2-Beg (1) or End (2), 3-Traj Length, 4-Frame#, 5-X Pos, 6-Y pos, 7-Placeholder for Scrum,
%8-Placeholder for # of scrums in frame, 9-#of scrums in interacting
%distance

for m=1:nvgranules
    temp = granules(granules(:,4) == vgranules(m),:);
    potinter(2*(m-1)+1,1) = temp(1,4);
    potinter(2*(m-1)+1,2) = 1;
    potinter(2*(m-1)+1,3) = size(temp,1);
    potinter(2*(m-1)+1,4) = temp(1,3);
    potinter(2*(m-1)+1,5) = temp(1,1);
    potinter(2*(m-1)+1,6) = temp(1,2);
    
    potinter(2*(m-1)+2,1) = temp(1,4);
    potinter(2*(m-1)+2,2) = 2;
    potinter(2*(m-1)+2,3) = size(temp,1);
    potinter(2*(m-1)+2,4) = temp(end,3);
    potinter(2*(m-1)+2,5) = temp(end,1);
    potinter(2*(m-1)+2,6) = temp(end,2);
end

waitbar(0,wh,['Granules Found: ' num2str(nvgranules) ' Analyzing Scrums'])

vscrums = zeros(max(scrums(:,4)),2); 
vscrums(:,1) = 1:max(scrums(:,4));

for m=1:max(scrums(:,4))
    temp = scrums(scrums(:,4) == m,:);
    if size(temp,1) >= validitylength
        vscrums(m,2) = 1;
    end
end

vscrums = vscrums(vscrums(:,2) == 1,1);
nvscrums = length(vscrums);

waitbar(0,wh,{['Granules Found: ' num2str(nvgranules) ' Scrums Found: ' num2str(nvscrums)]; ...
    ['Analyzing Potential Interaction: ' num2str(0) ' of ' num2str(nvgranules*2)]})

update = 100;
for m=1:2*nvgranules
    
    %find the scrums in the frame of interaction
    tscrum = scrums(scrums(:,3) == potinter(m,4),:);
    ntscrum = size(tscrum,1);
    potinter(m,8) = ntscrum; 

    %now test each scrum for interaction
    dis = zeros(ntscrum,3);
    dis(:,1) = tscrum(:,4);
    dis(:,2) = sqrt((tscrum(:,1)-potinter(m,5)).^2 + (tscrum(:,2)-potinter(m,6)).^2);
    dis(:,3) = dis(:,2) <= radius; 
    
    filtdist = dis(dis(:,3) == 1,:);
    
    if size(filtdist,1) == 1 
        potinter(m,7) = filtdist(1,1);
        potinter(m,9) = size(filtdist,1);
    elseif size(filtdist,1) > 1
        disp(['Error in interaction ' num2str(m) '. ' ...
            num2str(size(filtdist,1)) ' scrums are in radius. Recording minimum distance.']);
        [~,idx] = min(filtdist(:,2));
        potinter(m,7) = filtdist(idx,1);
        potinter(m,9) = size(filtdist,1);
    end        
    
    if m == update
        waitbar(m/(2*nvgranules),wh,{['Granules Found: ' num2str(nvgranules) ' Scrums Found: ' num2str(nvscrums)]; ...
            ['Analyzing Potential Interaction: ' num2str(m) ' of ' num2str(nvgranules*2)]})
        update = update+100;
    end
    
end

interactions = potinter(potinter(:,7) > 0,:);

save(['Working/' namestr '_radius_' num2str(radius) '_granulescruminter.mat'], ...
    'validitylength','radius','interactions','potinter','vscrums','scrums');

close(wh);