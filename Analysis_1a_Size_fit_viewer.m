function Analysis_1a_Size_fit_viewer()
global gridX gridY ptImg results boxsize t tfit pfit
close all



[filename, pathname] = uigetfile({'*.mat','.mat files'},'select .mat file');

if isequal(filename,0)
    return;
end

load([pathname filename],'ptImg','results','boxsize','pixelsize');

[gridX,gridY] = meshgrid(1:0.25:boxsize+1);
pfit = 1:0.25:boxsize+1;
t = 0:pixelsize:boxsize*pixelsize;
tfit = 0:pixelsize/4:boxsize*pixelsize;
h = figure(1);

display_figure(h,1);

answer = inputdlg(['Point Number (1-' num2str(size(results,1))  '):'],'Select Point Number');

while(~isempty(answer))
    input = round(str2double(answer{1}));
    if isnan(input) || input < 1 || input > size(results,1)
        eh = errordlg('Input a valid number.');
        uiwait(eh);
    else
        display_figure(h,input);
    end
    
    answer = inputdlg(['Point Number (1-' num2str(size(results,1))  '):'],'Select Point Number');    
end

close all

end

function display_figure(h,m)
global gridX gridY ptImg results boxsize t tfit pfit
set(0,'CurrentFigure',h);
result = results(m,2:7);
fres = result(1).*exp(- ( (gridX-result(2)).^2./(2*result(4)^2) + (gridY-result(3)).^2./(2*result(5)^2) ) ) + result(6);

xprofile = ptImg(round2(result(3)),1:boxsize+1,m);
yprofile = ptImg(1:boxsize+1,round2(result(2)),m);
xfit = result(1).*exp(- ( ((pfit)-result(2)).^2./(2*result(4)^2) ) ) + result(6);
yfit = result(1).*exp(- ( ((pfit)-result(3)).^2./(2*result(5)^2) ) ) + result(6);
I = double(ptImg(:,:,m));
I = I-min(min(I));
I = I./max(max(I));
fres = fres-min(min(fres));
fres = fres./max(max(fres));
subplot(2,2,1)
imshow(I,'Colormap',parula);
subplot(2,2,2)
imshow(fres,'Colormap',parula)
%contourf(flipud(fres));
subplot(2,2,3)
plot(t,xprofile,'o',tfit,xfit)
xlabel('position (nm)')
ylabel('intensity')
title(['X profile point: ' num2str(m)])
%set(gca,'XAxisLocation','top')
subplot(2,2,4)
plot(t,yprofile,'o',tfit,yfit)
xlabel('position (nm)')
ylabel('intensity')
title(['Y profile point: ' num2str(m)])
%set(gca,'YAxisLocation','right')
end

function t = round2(x)
t = round(x);

if(t<1)
    t=1;
end

end