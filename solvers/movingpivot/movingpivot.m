function [dxdt] = movingpivot(t, x, PD)
%% [dxdt] = movingpivot(t, x, PD) Moving Pivot method for Nucleation and Growth
% Solves the PBE according to the Moving Pivot method (cf. e.g.  Kumar, S.; Ramkrishna, D. Chemical Engineering Science 1997, 52, 4659–4679.)
% Needs to solve ODE's for number of particles (N), pivot length (y), boundaries (boundaries) and concentration,
% i.e. the number of ODE's to be solved is 3*ngrid-2+1 and has an event listener.
% In general, the moving pivot method yields the most accurate results but
% has the highest computational cost. The method is particularly well
% suited for problems with discontinuous distributions and geometric grids.
% If you are unhappy with the result consider the following options:
% - Increase the number of grid points
% - Decrease reltol {1e-6} and abstol {1e-6} [sol_options]
% - Reduce dL (criterion for bin addition) {10} [sol_options]
% - Use another method

nBins = (length(x)-2)/3; % number of bins   

% Current concentration
c = x(end); % 
    
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
y = x(nBins+1:2*nBins); %pivot sizes
boundaries = x(2*nBins+1:3*nBins+1); %boundaries
Dy = diff(boundaries);Dy = Dy(:);

if nargin(PD.nucleationrate)>3
    dist = Distribution(y,N./Dy,boundaries);
    J = PD.nucleationrate(S,T,t,dist);
else
    J = PD.nucleationrate(S,T,t);
end

Gy = PD.growthrate(S,T,y); % growth rate for pivots
Gboundaries = PD.growthrate(S,T,boundaries); % growth rate for boundaries

Gboundaries(boundaries<=0) = 0;

dNdt = [J; zeros(nBins-1,1)]-N(1:nBins)/m*Q; % change in number (per mass medium): nucleation - dilution
dcdt = -3*PD.rhoc*PD.kv*sum(y.^2.*Gy.*N)-c/m*Q-J*y(1)^3*PD.kv*PD.rhoc;

dxdt = [dNdt; Gy; Gboundaries; dcdt;];


end