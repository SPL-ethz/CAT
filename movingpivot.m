function [dxdt] = movingpivot(~, x, PD)
% movingpivot method for nucleation and growth, the 

    nBins = (length(x)-4)/3; % number of bins   
    
    N = x(1:nBins); %particle numbers
    p = x(nBins+1:2*nBins); %pivot sizes
    L = x(2*nBins+1:3*nBins+1); %boundaries
    c = x(3*nBins+2); %solution concentration 
    T = x(3*nBins+3); %temperature
    V = x(3*nBins+4); %reactor volume (non-constant in case of AS crystallization)
    
    J = PD.nucleationrate(c,T);
    G = PD.growthrate(c,T,[p; L]);
    
    dNdt = [J; -N(2:nBins)/V*PD.ASadditionrate];
    dpdt = [0.5*G(1); G(2:nBins)'];
    dLdt = [0; G(nBins+2:2*nBins+1)'];
    
    dcdt = -3*PD.rhoc*PD.kv*sum(p.^2.*dpdt.*N)-c/V*PD.ASadditionrate;
    dTdt = -PD.coolingrate;
    dVdt = PD.ASadditionrate;
    dxdt = [dNdt; dpdt; dLdt; dcdt; dTdt; dVdt];
end