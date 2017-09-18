clearvars;

namestr = '170427_3B11M_P13_Plate2a_Top'; 


load(['Working/' namestr '_sizeinter.mat']); 
load(['Working/' namestr '_msd_new.mat']); 

msd_nan = msd;
msd_nan(msd_nan == 0) = NaN;


T = table((1:length(ids))',ids(1,:)',ids(2,:)',ids(3,:)','VariableNames',{'rownum','id','sizeclass','length'});

csvwrite(['Output/' namestr '_msdfull_R.csv'],msd_nan);
writetable(T,['Output/' namestr '_msdids_R.csv']);