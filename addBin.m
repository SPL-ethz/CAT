function [x] = addBin(x)
% Adds a new bin to the distribution, so that we can continue integrating
% without producing too much of an error in the mass balance.
    nBins = (length(x)-4)/3;
    
    N = x(1:nBins); %particle numbers
    p = x(nBins+1:2*nBins); %pivot sizes
    L = x(2*nBins+1:3*nBins+1); %boundaries
    c = x(3*nBins+2); %solution concentration
    T = x(3*nBins+3); %solution concentration
    V = x(3*nBins+4); %solution concentration
    
    %shift old contents and add new bin
    N = [0; N];
    p = [0; p];
    L = [0; L];    
    
    %re-assemble x
    x = [N; p; L; c; T; V];
end