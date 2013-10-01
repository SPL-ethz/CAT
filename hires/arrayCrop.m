function y = arrayCrop(x,cropSize)

sfx = size(x);sfx(sfx==1) = [];
nDim = length(sfx);

ind = {};

for k = 1:nDim
   ind = [ind cropSize(1,k):cropSize(2,k)];
end

y = x(ind{:});
