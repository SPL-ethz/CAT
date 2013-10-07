function [GI] = boundingBoxFinder(fstar,mlim)

fx1 = cumsum(abs(fstar));

a0 = find(fx1./fx1(end)>mlim,1,'first');

a0(isempty(a0)) = 1;

af = find(1-fx1./fx1(end)<mlim,1,'first');
af(isempty(af)) = length(fstar(3:end-1,1,1));

GI(1,1) = min([a0 af-1]);
GI(1,2) = af;
