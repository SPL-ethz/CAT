function [mbflag] = solver_centraldifference(O)
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

if ~isempty(O.sol_options) && ~isempty(find(strcmpi(O.sol_options,'massbalTol'),1))
    massbalTol = O.sol_options(find(strcmpi(O.sol_options,'massbalTol'),1)+1);
else
    massbalTol = 0.05; % massbalance error tolerance
end

if isempty(O.sol_options) || isempty(O.sol_options{1})
    options = odeset('Events',@(t,x) Event(t,x,massbalTol,O),'reltol',1e-6);
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

if O.calc_time(end)<O.sol_time(end) 
    mbflag = 1;
else
    mbflag = 0;
end


end % function

function dXdt = centraldifference_ode(t,X,O)

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
S = c/evalanonfunc(O.solubility,T,xm);

% Current mass flow rate antisolvent (evaluated using simplistic FFD)
Q = (O.ASprofile(t+1e-6)-O.ASprofile(t))/1e-6;
if isnan(Q)
    Q = (O.ASprofile(t)-O.ASprofile(t-1e-6))/1e-6;
end

%% Growth rate evaluation
G = evalanonfunc(O.growthrate, S, T, y(:), t );

Ga = evalanonfunc(O.growthrate, S, T, ya, t );
Fa = F( [2:end end] );

Gb = evalanonfunc(O.growthrate, S, T, yb, t );
Fb = F( [1 1:end-2 end-2] );

% nucleation
if nargin(O.nucleationrate)>2
    dist = Distribution(y,F,O.init_dist.boundaries);
else
    dist = [];
end

J = evalanonfunc(O.nucleationrate, S, T, dist, t );

% nucleation
Fb(1) = J/G(1);

% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';
dF(isnan(dF)) = 0;

% concentration
dc = -3*O.kv*O.rhoc*sum(G(:).*F(:).*Dy(:).*y(:).^2)-c/m*Q-J*y(1)^3*O.kv*O.rhoc;

dXdt = [dF(:)-Q*F(:)/m; dc];

% if GUI is activated, update progress bar
if findall(0,'name','Looking at CATs')
            
    if isempty(O.tNodes)
        tFinal = O.sol_time(end);
    else
        tFinal = O.tNodes(end);
    end
    if floor(t/tFinal/0.05)>str2num(get(gca,'tag'))
        fill([0 t/tFinal t/tFinal 0],[0 0 1 1],'c','edgecolor','none')
        delete(findall(gcf,'type','text'))
        text(0.44,0.5,[num2str(floor(t/tFinal*100),'%2d'),'%'])
        set(gca,'tag',num2str(floor(t/tFinal/0.05)))
        drawnow
    end
end


end % function

%% Function EventBin

function [value,isterminal,direction] = Event(t,x,massbalTol,O)

m = O.init_massmedium+(O.ASprofile(t)-O.ASprofile(0));

value(1) = massbalTol-...
    abs(((O.init_conc*O.init_massmedium+moments(O.init_dist,3)*O.kv*O.rhoc*O.init_massmedium)-(x(end)*m+sum(x(1:end-1).*diff(O.init_dist.boundaries(:)).*O.init_dist.y(:).^3)*O.kv*O.rhoc*m))/...
    (O.init_conc*O.init_massmedium+moments(O.init_dist,3)*O.kv*O.rhoc*O.init_massmedium)); % current massbalance error

value = value(:);
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end