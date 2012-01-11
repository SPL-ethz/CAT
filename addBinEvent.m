function [value,isterminal,direction] = addBinEvent(~,x,dL)
% Locate the time when the cell where nucleation happens becomes too big.

    nBins = (length(x)-2)/3;
    
    L = x(2*nBins+1:2*nBins+2);
        
    value = dL - L(2)+L(1);     % Detect when nucleation bin becomes too broad
    isterminal = 1;   % Stop the integration
    direction = -1;   % Negative direction only
end