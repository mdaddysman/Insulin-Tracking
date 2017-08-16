function data = histogram2pts(X,varargin)
%HISTOGRAMPTS Summary of this function goes here
%   Detailed explanation goes here

figure(100)
h = histogram(X,varargin{:});
be = h.BinEdges;
bc = h.Values; 
bw = h.BinWidth; 
bcen = be(1:end-1)+bw/2;
close 100

data = [bcen',bc'];

end

