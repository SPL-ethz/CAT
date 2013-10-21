function [GI] = boundingBoxFinder(F,mlim)

Fcum = cumsum(abs(F));

a0 = find(Fcum./Fcum(end)>mlim,1,'first');
a0(isempty(a0)) = 1; 

af = find(1-Fcum./Fcum(end)<mlim,1,'first');
af(isempty(af)) = length(F(3:end-1,1,1));

GI(1,1) = min([a0 max([1 af-1])]);
GI(1,2) = max([2 af]);
