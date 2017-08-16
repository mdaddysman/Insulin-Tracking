clearvars;
namestr = '170721imaris'; 

tracked = csvread([namestr '.csv']); 

%this script adds 5 columns & combines the tracked variables
%col 6 is a small (1) or large (2) particle trajectory -> all small 
%col 7 is the distance to the closest particle in pixel units
%col 8 is the size of the particle from a gaussian fit 
%col 9 is a boolean if the size is "valid" (near the center of the image
%box)
%col 10 is the ID for looking up the fit results in the viewer (1a) 
%Position X	Position Y	Position Z	Unit	Category	Collection	Time	TrackID	ID

tracked(:,5) = tracked(:,5) - 10^9 + 1;



%expand the size of the matrix
tracked(:,6:10) = zeros([size(tracked,1) 5]);
tracked(:,6) = ones([size(tracked,1) 1]);



%restore order - sort by traj, frame, large small 
tracked = sortrows(tracked,[5 4]);
save([namestr '_sizeinter.mat'],'tracked');
