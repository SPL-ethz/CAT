function [xout] = addBin(xin)
% Adds a new bin to the distribution, so that we can continue integrating
% without producing too much of an error in the mass balance.

nBins = (length(xin)-2)/3; % number of bins in distribution

N = xin(1:nBins); %particle numbers
y = xin(nBins+1:2*nBins); % pivot sizes
boundaries = xin(2*nBins+1:3*nBins+1); %boundaries
c = xin(3*nBins+2); % concentration

%insert new, empty bin at beginning (we only support nucleation)
N = [0; N];
y = [0; y];
boundaries = [0; boundaries];    

% assemble new x
xout = [N; y; boundaries; c;];
end