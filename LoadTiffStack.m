function I = LoadTiffStack(filename)
%LOADTIFFSTACK Load a multipage tiff stack into memory as a 3D vector

info = imfinfo(filename);

width = info(1).Width;
height = info(1).Height;
bitdepth = info(1).BitDepth;
images = size(info,1);

switch bitdepth
    case 8
        I = zeros([height width images],'uint8');
    case 16
        I = zeros([height width images],'uint16');
    otherwise
        I = zeros([height width images],'double');
end

for m=1:size(info,1)
    I(:,:,m) = imread(filename,'Info',info,'Index',m);
end

end

