function y = gauss1D(x,t)
%GAUSS1D Summary of this function goes here
%   Detailed explanation goes here

A = x(1);
xc = x(2);
xw = x(3);
C = x(4);

y = A.*exp(- ( (t-xc).^2./(2*xw^2) ) ) + C;

end

