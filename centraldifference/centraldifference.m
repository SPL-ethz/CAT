function dXdt = centraldifference(t,X,PD)

% Current solvent + antisolvent mass
m = PD.init_massmedium+(PD.ASprofile(t)-PD.ASprofile(0));

% Current mass fraction antisolvent
xm = PD.ASprofile(t)/m;

% Grid
y = PD.init_dist.y;
Dy = diff([0 y(:)']);
ya = y([2:end end]);
yb = y([1 1:end-2 end-2]);

% Current concentration and distribution
F = X(1:length(y))'; % distribution
c = X(end); % concentration

% Current temperature
T = PD.Tprofile(t);

% Current supersaturation
S = c/PD.solubility(T,xm);

% Current mass flow rate antisolvent (evaluated using simplistic FD)
Q = (PD.ASprofile(t+1e-12)-PD.ASprofile(t-1e-12))/2e-12;

G = PD.growthrate(S,T,y(:));

Ga = PD.growthrate(S, T, ya );
Fa = F( [2:end end] );

Gb = PD.growthrate(S, T, yb );
Fb = F( [1 1:end-2 end-2] );

% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';
dc = -3*PD.kv*PD.rhoc*sum(G.*F(:).*Dy(:).*y(:).^2)-c/m*Q;

% nucleation
J = PD.nucleationrate(S,T);
X(1) = J/G(1);

dXdt = [dF-Q*X(1:end-1)/m; dc];

end % function