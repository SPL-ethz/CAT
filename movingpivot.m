function [dxdt] = movingpivot(~, x, PD)
% some commments 
    nBins = (length(x)-2)/3;
    T = PD.init_temp;
    
    N = x(1:nBins); %particle numbers
    p = x(nBins+1:2*nBins); %pivot sizes
    L = x(2*nBins+1:3*nBins+1); %boundaries
    c = x(3*nBins+2); %solution concentration  
    
    J = PD.nucleationrate(c,T);
    G = PD.growthrate(c,T,[p; L]);
    
    dNdt = [J; zeros(nBins-1,1)];
    dpdt = [0.5*G(1); G(2:nBins)'];
    dLdt = [0; G(nBins+2:2*nBins+1)'];
    
    dcdt = -3*PD.rhoc*PD.kv*sum(p.^2.*dpdt.*N);
    dxdt = [dNdt; dpdt; dLdt; dcdt];
end