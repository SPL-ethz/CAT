function dXdt = centraldifference(t,X,PD)

% Current solvent + antisolvent mass
m = PD.init_massmedium+(PD.ASprofile(t)-PD.ASprofile(0));

% Current mass fraction antisolvent
xm = PD.ASprofile(t)/m;

% Current concentration
c = X(end);

% Grid
y = PD.init_dist.y;
ya = y([2:end end]);
yb = y([1 1:end-2 end-2]);

% Current temperature
T = PD.Tprofile(t);

% Current supersaturation
S = c/PD.solubility(T,xm)

% Current mass flow rate antisolvent (evaluated using simplistic FD)
Q = (PD.ASprofile(t+1e-12)-PD.ASprofile(t-1e-12))/2e-12;
    

F = X(1:length(y))';

Ga = PD.growthrate(S, T, ya );
Fa = F( [2:end end] );

Gb = PD.growthrate(S, T, yb );
Fb = F( [1 1:end-2 end-2] );

% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';
dc = -3*PD.kv*PD.rhoc*sum(PD.growthrate(S,T,y(:)).*X(1:end-1).*y(:).^2)-c/m*Q;

dXdt = [dF-Q*X(1:end-1)/m; dc];

end % function