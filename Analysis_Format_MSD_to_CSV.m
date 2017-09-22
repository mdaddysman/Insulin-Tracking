clearvars;

namestr = '161101_3D3_P11_2d'; 


load(['Working/' namestr '_sizeinter.mat']); 
load(['Working/' namestr '_msd_new.mat']); 

msd_nan = msd;
msd_nan(msd_nan == 0) = NaN;


T = table((1:length(ids))',ids(1,:)',ids(2,:)',ids(3,:)','VariableNames',{'rownum','id','sizeclass','length'});
Tpos = table(tracked(:,3),tracked(:,1),tracked(:,2),tracked(:,4),tracked(:,5), ...
    'VariableNames',{'time','x','y','id','class'});

csvwrite(['Output/' namestr '_msdfull_R.csv'],msd_nan);
writetable(T,['Output/' namestr '_msdids_R.csv']);
writetable(Tpos,['Output/' namestr '_trajpos_R.csv']);