function F = gauss2D(x)
%GAUSS2D Summary of this function goes here
%   Detailed explanation goes here
global gridX gridY fitpicture 

A = x(1);
xc = x(2);
yc = x(3);
xw = x(4);
yw = x(5);
C = x(6);

y = A.*exp(- ( (gridX-xc).^2./(2*xw^2) + (gridY-yc).^2./(2*yw^2) ) ) + C;

F = y - fitpicture; 

end

