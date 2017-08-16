function corrected_p = SPIFF_sort(p_in,SPIFF_x,SPIFF_y)

%%%%%%%%%%%%SPIFF correction %%%%%%%%%%%%%%%%%%
%%%Step 1 - check SPIFF data

%%%% R = 2   %%%%%
inds = find(p_in(:,1)>0);


%Now calculate the positive and negative probability density distributions
%(P(x_e) in the manuscript - eq. (7));
SPIFF_x_plus = SPIFF_x(SPIFF_x>=0.5);
SPIFF_y_plus = SPIFF_y(SPIFF_y>=0.5);
SPIFF_x_minus = SPIFF_x(SPIFF_x<0.5);
SPIFF_y_minus = SPIFF_y(SPIFF_y<0.5);

p_x_plus = sort(SPIFF_x_plus);
p_x_minus = sort(SPIFF_x_minus);
p_y_plus = sort(SPIFF_y_plus);
p_y_minus = sort(SPIFF_y_minus);

%The SPIFF corrected positions go here:
corrected_p = zeros(size(p_in));
  
%Seperates the positive and negative parts of the SPIFF in both x and y.
% for each one of the biased points, calculates the amount of points
% below its value.


for ww = 1:length(p_in)
    %Treat position list for particle 1
    if p_in(ww,1) == 0
        continue
    end
    x1_curr = p_in(ww,1);
    y1_curr = p_in(ww,2);
    SPIFF_x_curr = x1_curr - floor(x1_curr);
    SPIFF_y_curr = y1_curr - floor(y1_curr);  
    if SPIFF_x_curr>=0.5
        %a. Get percentile value:
        tmp_x = sum(p_x_plus<=SPIFF_x_curr)/length(p_x_plus);
        %b. Bring back to unmodolu'ed position:
        new_SPIFF_x = tmp_x/2; %if the position is positive (above the floor) it stays.
        corrected_p(ww,1) = ceil(x1_curr) + new_SPIFF_x;
    else
         %a. Get percentile value (how many particles are closer to center):
        tmp_x = sum((p_x_minus)>=(SPIFF_x_curr))/length(p_x_minus);
        %b. Bring back to unmodolu'ed position:
        new_SPIFF_x = tmp_x/2;   %If the position is negative (below the floor) it is moved by 1. 
        corrected_p(ww,1) = ceil(x1_curr) - new_SPIFF_x;
    end
    if SPIFF_y_curr>=0.5
        %a. Get percentile value:
        tmp_y = sum(p_y_plus<=SPIFF_y_curr)/length(p_y_plus);
        %b. Bring back to unmodolu'ed position:
        new_SPIFF_y = tmp_y/2; %if the position is positive (above the floor) it stays.
        corrected_p(ww,2) = ceil(y1_curr) + new_SPIFF_y;
    else
         %a. Get percentile value (how many particles are closer to center):
        tmp_y = sum((p_y_minus)>=(SPIFF_y_curr))/length(p_y_minus);
        %b. Bring back to unmodolu'ed position:
        new_SPIFF_y = tmp_y/2;   %If the position is negative (below the floor) it is moved by 1. 
        corrected_p(ww,2) = ceil(y1_curr) - new_SPIFF_y;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%