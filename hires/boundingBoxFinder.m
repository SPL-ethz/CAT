function [GI] = boundingBoxFinder(fstar,ndim,mlim,volOn,V)
% keyboard
if nargin>3 && volOn == 1
    V = V;
else
    V = 1;
end

if ndim == 1
    fx1 = cumsum(abs(fstar).*V);
%     fx1(1:2) = [];fx2(1:2) = [];
%     fx1(end) = [];fx2(end) = [];
%                 tic
    a0 = find(fx1./fx1(end)>mlim,1,'first');

    a0(isempty(a0)) = 1;

    af = find(1-fx1./fx1(end)<mlim,1,'first');
    af(isempty(af)) = length(fstar(3:end-1,1,1));

    GI(1,1) = min([a0 af-1]);
    GI(1,2) = af;

elseif ndim == 2
    fx1 = sum(cumsum(abs(fstar).*V),2);
    fx2 = sum(cumsum(abs(fstar).*V,2));

%     fx1(1:2) = [];fx2(1:2) = [];
%     fx1(end) = [];fx2(end) = [];
%                 tic
    a0 = find(fx1./fx1(end)>mlim,1,'first');

    a0(isempty(a0)) = 1;

    af = find(1-fx1./fx1(end)<mlim,1,'first');
    af(isempty(af)) = length(fstar(3:end-1,1,1));

    GI(1,1) = min([a0 af-1]);
    GI(1,2) = af;
%                 GI(1,2) = min([af length(G.dim1)]);
    b0 = find(fx2./fx2(end)>mlim,1,'first');
    b0(isempty(b0)) = 1;
    bf = find(1-fx2./fx2(end)<mlim,1,'first');
    bf(isempty(bf)) = length(fstar(1,3:end-1,1));

    GI(2,1) = min([b0 bf-1]);
    GI(2,2) = bf;

elseif ndim == 3
%                 keyboard
    fx1 = sum(sum(cumsum(abs(fstar).*V  ),2),3);
    fx2 = sum(sum(cumsum(abs(fstar).*V,2),3));
    fx3 = sum(sum(cumsum(abs(fstar).*V,3)));

%     fx1(1:2) = [];fx2(1:2) = [];fx3(1:2) = [];
%     fx1(end) = [];fx2(end) = [];fx3(end) = [];

    a0 = find(fx1./fx1(end)>mlim,1,'first');
    a0(isempty(a0)) = 1;

%                 keyboard
    af = find(1-fx1./fx1(end)<mlim,1,'first');
    af(isempty(af)) = length(fstar(3:end-1,1,1));

    GI(1,1) = a0;
    GI(1,2) = af;


    b0 = find(fx2./fx2(end)>mlim,1,'first');
    b0(isempty(b0)) = 1;
    bf = find(1-fx2./fx2(end)<mlim,1,'first');
    bf(isempty(bf)) = length(fstar(1,3:end-1,1));

    GI(2,1) = b0;
    GI(2,2) = bf;

    c0 = find(fx3./fx3(end)>mlim,1,'first');
    c0(isempty(c0)) = 1;
    cf = find(1-fx3./fx3(end)<mlim,1,'first');
    cf(isempty(cf)) = length(fstar(1,1,3:end-1));

    GI(3,1) = c0;
    GI(3,2) = cf;
end