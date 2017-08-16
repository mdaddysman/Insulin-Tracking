function DisplayStack(I,varagin)
%DISPLAYRGBSTACK Makes a figure that has a slider bar to change a stack of
%RGB images
switch nargin
    case 2
        if ischar(varagin) == 1
            fh = figure('Name',varagin);
        else
            fh = figure;
        end
    case 1
        fh = figure;
    otherwise
        disp('Only use 1 or 2 inputs.');
        return;
end

switch length(size(I))
    case 3

        ih = imshow(I(:,:,1));
        
        pos = get(fh,'Position');
        set(fh,'Resize','off');
        
        nth = uicontrol('Style','text','String','1','Position',[5 5 45 30],'FontSize',24,'FontWeight','bold');
        majorstep = 5/(size(I,3) - 1);
        minorstep = 1/(size(I,3) - 1);
        uicontrol('Style','slider','Min',1,'Max',size(I,3),'Value',1,'Position',[55 5 pos(3)-65 30], ...
            'SliderStep',[minorstep majorstep],'Callback',{@GrayscaleStack_SliderControl,nth,ih,I});
    case 4
        ih = imshow(I(:,:,:,1));
        
        pos = get(fh,'Position');
        set(fh,'Resize','off');
        
        nth = uicontrol('Style','text','String','1','Position',[5 5 45 30],'FontSize',24,'FontWeight','bold');
        majorstep = 5/(size(I,4) - 1);
        minorstep = 1/(size(I,4) - 1);
        uicontrol('Style','slider','Min',1,'Max',size(I,4),'Value',1,'Position',[55 5 pos(3)-65 30], ...
            'SliderStep',[minorstep majorstep],'Callback',{@RGBStack_SliderControl,nth,ih,I});
    otherwise
        disp('Image stack must be 3 (grayscale) or 4 (RGB) dimensional.');
        close(fh);
        return;
end

end

function RGBStack_SliderControl(src,~,nth,ih,I)
value = get(src,'Value');
value = round(value);

set(nth,'String',num2str(value));
set(src,'Value',value);

set(ih,'CData',I(:,:,:,value));

end

function GrayscaleStack_SliderControl(src,~,nth,ih,I)
value = get(src,'Value');
value = round(value);

set(nth,'String',num2str(value));
set(src,'Value',value);

set(ih,'CData',I(:,:,value));

end