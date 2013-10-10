function [xout] = removeBin(xin)
% Removes a bin from the distribution when necessary (dissolution)

nBins = (length(xin)-2)/3; % number of bins in distribution

N = xin(1:nBins); %particle numbers
y = xin(nBins+1:2*nBins); % pivot sizes
boundaries = xin(2*nBins+1:3*nBins+1); %boundaries
c = xin(3*nBins+2); % concentration

% remove first bin
N = N(2:end);
y = y(2:end);
boundaries = boundaries(2:end);    

% assemble new x
xout = [N; y; boundaries; c;];
end