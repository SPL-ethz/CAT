function [xout] = addBin(xin)
% Adds a new bin to the distribution, so that we can continue integrating
% without producing too much of an error in the mass balance.
    nBins = (length(xin)-2)/3;
    
    N = xin(1:nBins); %particle numbers
    p = xin(nBins+1:2*nBins); %pivot sizes
    L = xin(2*nBins+1:3*nBins+1); %boundaries
    c = xin(3*nBins+2); %solution concentration
    
    %shift old contents and add new bin
    N = [0; N];
    p = [0; p];
    L = [0; L];    
    
    %re-assemble x
    xout = [N; p; L; c;];
end