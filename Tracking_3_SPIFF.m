namestr = '170427_3B11M_P13_Plate2a_Top';

load([namestr '_filtered.mat']);

data_s = smalldata(:,2:3);
mp_totx = data_s(:,1)-floor(data_s(:,1));
mp_toty = data_s(:,2)-floor(data_s(:,2));
SPIFF_s = SPIFF_sort(data_s,mp_totx,mp_toty);
smalldata(:,2:3) = SPIFF_s;


data_l = largedata(:,2:3);
mp_totx = data_l(:,1)-floor(data_l(:,1));
mp_toty = data_l(:,2)-floor(data_l(:,2));
SPIFF_l = SPIFF_sort(data_l,mp_totx,mp_toty);
largedata(:,2:3) = SPIFF_l;

save([namestr '_SPIFF.mat'],'smalldata','largedata');

