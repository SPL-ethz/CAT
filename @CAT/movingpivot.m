function [dxdt] = movingpivot(O, t, x)
%% [dxdt] = movingpivot(t, x, O) Moving Pivot method for Nucleation and Growth
% Solves the PBE according to the Moving Pivot method (cf. e.g.  Kumar, S.; Ramkrishna, D. Chemical Engineering Science 1997, 52, 4659�4679.)
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
m = O.init_massmedium+(O.ASprofile(t)-O.ASprofile(0));

% Current mass fraction antisolvent
xm = O.ASprofile(t)/m;

% Current temperature
T = O.Tprofile(t);

% Current supersaturation
S = c/O.solubility(T,xm);

% Current mass flow rate antisolvent (evaluated using simplistic FD)
Q = (O.ASprofile(t+1e-6)-O.ASprofile(t))/1e-6;
if isnan(Q)
    Q = (O.ASprofile(t)-O.ASprofile(t-1e-6))/1e-6;
end


N = x(1:nBins); N = N(:); %particle numbers
y = x(nBins+1:2*nBins); %pivot sizes
boundaries = x(2*nBins+1:3*nBins+1); %boundaries
Dy = diff(boundaries);Dy = Dy(:);

if nargin(O.nucleationrate)>3
    dist = Distribution(y,N./Dy,boundaries);
    J = O.nucleationrate(S,T,t,dist);
else
    J = O.nucleationrate(S,T,t);
end

Gy = O.growthrate(S,T,y); % growth rate for pivots
Gboundaries = O.growthrate(S,T,boundaries); % growth rate for boundaries

Gboundaries(boundaries<=0) = 0;

dNdt = [J; zeros(nBins-1,1)]-N(1:nBins)/m*Q; % change in number (per mass medium): nucleation - dilution
dcdt = -3*O.rhoc*O.kv*sum(y.^2.*Gy.*N)-c/m*Q-J*y(1)^3*O.kv*O.rhoc;

dxdt = [dNdt; Gy; Gboundaries; dcdt;];


end

%% Function addBin

function [xout] = addBin(xin)
% Adds a new bin to the distribution, so that we can continue integrating
% without producing too much of an error in the mass balance.

nBins = (length(xin)-2)/3; % number of bins in distribution

N = xin(1:nBins); %particle numbers
y = xin(nBins+1:2*nBins); % pivot sizes
boundaries = xin(2*nBins+1:3*nBins+1); %boundaries
c = xin(3*nBins+2); % concentration

%insert new, empty bin at beginning (we only support nucleation)
N = [0; N];
y = [0; y];
boundaries = [0; boundaries];

% assemble new x
xout = [N; y; boundaries; c;];
end

%% Function EventBin

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

%% Ffunction removeBin

function [xout] = removeBin(xin)
% Removes a bin from the distribution when necessary (dissolution)

nBins = (length(xin)-2)/3; % number of bins in distribution

N = xin(1:nBins); %particle numbers
y = xin(nBins+1:2*nBins); % pivot sizes
boundaries = xin(2*nBins+1:3*nBins+1); %boundaries
c = xin(3*nBins+2); % concentration

% remove first bin
N = N(2:end);
y = y(2:end);
boundaries = boundaries(2:end);

% assemble new x
xout = [N; y; boundaries; c;];
end

