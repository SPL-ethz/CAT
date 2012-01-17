%% solubility as function of temperature

% cs: solubility in kg solute per kg solvent [kg/m^3]
% T: Temperature [K]

function [cs] = SolubilityAlphaLGLU(T)
    %solubility data from:
    % Y. Mo et al., Fluid Phase Equilibria 300 (2011) 105-109
    x =  exp(-26.9e3/(8.3145.*T)+35.94/(8.3145)); %x is mole fraction
    M1 = 147.1;
    M2 = 18;
    xM = x*M1/(x*M1+(1-x)*M2);%xM is mass fraction
    cs = xM/(1-xM)*1000;    
end  