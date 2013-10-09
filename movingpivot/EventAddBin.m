function [value,isterminal,direction] = EventAddBin(~,x,dL)
% Event function 

nBins = (length(x)-2)/3;

boundaries = x(2*nBins+1:2*nBins+2);

value = dL - boundaries(1); % Detect when first bin becomes too broad (value <= 0)
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only
end