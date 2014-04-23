function solver_centraldifference(O)
%% solver_centraldifference
% Central Difference Method for Nucleation and Growth
% Solves the PBE according to a central differences method (cf. Wikipedia).
% Needs to solve ODE's for number of particles (N), pivot length (y), boundaries (boundaries) and concentration,
% i.e. the number of ODE's to be solved is ngrid-1+1.
% In general, this method yields reasonably accurate results with low
% computational cost.
% Central differences methods are good allrounders but are prone to
% oscillations when dealing with discontinuous distributions. Also, the
% mass balance error is considerably larger than with the moving pivot
% method.
% If you are unhappy with the result consider the following options:
% - Increase the number of grid points
% - Decrease reltol {1e-6} and abstol {1e-6} [sol_options]
% - Use another method


options = O.sol_options;
if isempty(O.sol_options)
    options = odeset(options,'reltol',1e-6);
end

X0 = [O.init_dist.F, O.init_conc];

solvefun = @(t,X) centraldifference_ode(t,X,O);

[SolutionTimes,X_out] = ode15s(solvefun , O.sol_time , X0 ,options);

% Transform result-arrays into appropriate output structure
SolutionConc = X_out(:,end);
SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  % Pre-Allocation for speed
for i = 1:length(SolutionTimes)
    SolutionDists(i) = Distribution( O.init_dist.y, X_out(i,1:length(O.init_dist.y)),O.init_dist.boundaries );
end % for

O.calc_time = SolutionTimes;
O.calc_dist = SolutionDists;
O.calc_conc = SolutionConc;

end % function

function centraldifference_ode(t,X,O)

% Current solvent + antisolvent mass
m = O.init_massmedium+(O.ASprofile(t)-O.ASprofile(0));

% Current mass fraction antisolvent
xm = O.ASprofile(t)/m;

% Grid
y = O.init_dist.y;
Dy = diff(O.init_dist.boundaries);
ya = y([2:end end]);
yb = y([1 1:end-2 end-2]);

% Current concentration and distribution
F = X(1:length(y))'; % distribution
c = X(end); % concentration

% Current temperature
T = O.Tprofile(t);

% Current supersaturation
S = c/O.solubility(T,xm);

% Current mass flow rate antisolvent (evaluated using simplistic FFD)
Q = (O.ASprofile(t+1e-6)-O.ASprofile(t))/1e-6;
if isnan(Q)
    Q = (O.ASprofile(t)-O.ASprofile(t-1e-6))/1e-6;
end

%% Growth rate evaluation
G = O.growthrate(S,T,y(:));

Ga = O.growthrate(S, T, ya );
Fa = F( [2:end end] );

Gb = O.growthrate(S, T, yb );
Fb = F( [1 1:end-2 end-2] );

% nucleation
if nargin(O.nucleationrate)>3
    dist = Distribution(y,F,O.init_dist.boundaries);
    J = O.nucleationrate(S,T,t,dist);
    % nucleation
    Fb(1) = J/G(1);
else
    J = O.nucleationrate(S,T,t);
    % nucleation
    Fb(1) = J/G(1);
end


% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';
dF(isnan(dF)) = 0;

% concentration
dc = -3*O.kv*O.rhoc*sum(G(:).*F(:).*Dy(:).*y(:).^2)-c/m*Q-J*y(1)^3*O.kv*O.rhoc;

dXdt = [dF(:)-Q*F(:)/m; dc];


end % function