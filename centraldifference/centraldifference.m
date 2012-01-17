function dXdt = centraldifference(~,X,PD)

c = X(end);
y = PD.init_dist.y;
T = PD.init_temp;
F = X(1:length(y))';

ya = y([2:end end]);
Ga = PD.growthrate(c, T, ya );
Fa = F( [2:end end] );

yb = y([1 1:end-2 end-2]);
Gb = PD.growthrate(c, T, yb );
Fb = F( [1 1:end-2 end-2] );

% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';

dXdt = [dF; 0];

end % function