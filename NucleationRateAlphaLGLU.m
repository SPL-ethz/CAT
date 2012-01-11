%%  in dependence of supersaturation, temperature and crystal size

% J: nucleation rate [#/(kg s)]
% T: Temperature [K]
% S: Supersaturation [-]

function [J] = NucleationRateAlphaLGLU(S,T)
    global parameters;
    %example data taken from Schöll et al., Chem. Eng. Technol., 2006, 29 (2),
    %257-264 (alpha L-glutamic acid)
    k = [1.3e27 163 1.4e8 10];
    J = k(1)*exp(-k(2)/log(S).^2) + k(3)*exp(-k(4)/log(S).^2);
end  