%% crystal growth rate in dependence of supersaturation, temperature and crystal size

% G: growth rate [m/s]
% T: Temperature [K]
% S: Supersaturation [-]

function [G] = GrowthRateAlphaLGLU(S,T,L)
    %example data taken from Schöll et al., Faraday Discuss., 2007, 136,
    %245-264 (alpha L-glutamic acid)

    k = [3.63e-4 3.72e3 5.42e4];
    G = k(1).*T.*exp(-k(2)./T)*(S-1).^(2/3).*log(S).^(1/6)*exp(-k(3)./(T.^2.*log(S)));
    G = G*ones(1,length(L)); %provide same growth rate for all the particle sizes
end  