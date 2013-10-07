function dXdt = centraldifference(t,X,PD)

% Current solvent + antisolvent mass
m = PD.init_massmedium+(PD.ASprofile(t)-PD.ASprofile(0));

% Current mass fraction antisolvent
xm = PD.ASprofile(t)/m;

% Grid
y = PD.init_dist.y;
% keyboard
Dy = diff(PD.init_dist.boundaries);
ya = y([2:end end]);
yb = y([1 1:end-2 end-2]);

% Current concentration and distribution
F = X(1:length(y))'; % distribution
c = X(end); % concentration

% Current temperature
T = PD.Tprofile(t);

% Current supersaturation
S = c/PD.solubility(T,xm);

% Current mass flow rate antisolvent (evaluated using simplistic FFD)
Q = (PD.ASprofile(t+1e-6)-PD.ASprofile(t))/1e-6;
% keyboard

G = PD.growthrate(S,T,y(:));

Ga = PD.growthrate(S, T, ya );
Fa = F( [2:end end] );

Gb = PD.growthrate(S, T, yb );
Fb = F( [1 1:end-2 end-2] );

% nucleation
if nargin(PD.nucleationrate)==3
    dist = Distribution(y,F,PD.init_dist.boundaries);
    J = PD.nucleationrate(S,T,dist);
    % nucleation
    Fb(1) = J/G(1);
else
    J = PD.nucleationrate(S,T);
    % nucleation
    Fb(1) = J/G(1);
end


% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';


% concentration
dc = -3*PD.kv*PD.rhoc*sum(G.*F(:).*Dy(:).*y(:).^2)-c/m*Q-J*y(1)^3*PD.kv*PD.rhoc;

dXdt = [dF(:)-Q*F(:)/m; dc];
% t/3600

end % function