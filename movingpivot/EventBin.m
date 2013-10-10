function [value,isterminal,direction] = EventBin(~,x,dL)
% Event function 
nBins = (length(x)-2)/3;

boundaries = x(2*nBins+1:2*nBins+2);

value(1) = dL - boundaries(1); % Detect when first bin becomes too broad (value <= 0)
value(2) = x(nBins+1);

value = value(:);
isterminal = ones(2,1);   % Stop the integration
direction = -ones(2,1);   % Negative direction only
end