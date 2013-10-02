function [dxdt] = movingpivot(t, x, PD)
% movingpivot method for nucleation and growth

nBins = (length(x)-2)/3; % number of bins   
c = x(end);
    
% Current solvent + antisolvent mass
m = PD.init_massmedium+(PD.ASprofile(t)-PD.ASprofile(0));

% Current mass fraction antisolvent
xm = PD.ASprofile(t)/m;

% Current temperature
T = PD.Tprofile(t);

% Current supersaturation
S = c/PD.solubility(T,xm);

% Current mass flow rate antisolvent (evaluated using simplistic FD)
Q = (PD.ASprofile(t+1e-6)-PD.ASprofile(t))/1e-6;
if isnan(Q)
    Q = (PD.ASprofile(t)-PD.ASprofile(t-1e-6))/1e-6;
end

    
N = x(1:nBins); N = N(:); %particle numbers
p = x(nBins+1:2*nBins); %pivot sizes
b = x(2*nBins+1:3*nBins+1); %boundaries

J = PD.nucleationrate(S,T);
Gp = PD.growthrate(S,T,p);
Gb = PD.growthrate(S,T,b);

dNdt = [J; zeros(nBins-1,1)]-N(1:nBins)/m*Q;    

dcdt = -3*PD.rhoc*PD.kv*sum(p.^2.*Gp.*N)-c/m*Q;
dxdt = [dNdt; Gp; Gb; dcdt;];

end