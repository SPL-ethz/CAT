function dXdt = centraldifference(t,X,y,growthrate)

F = X(1:length(y))';
c = X(end);

ya = y([2:end end]);
Ga = growthrate(c, ya );
Fa = F( [2:end end] );

yb = y([1 1:end-2 end-2]);
Gb = growthrate(c, yb );
Fb = F( [1 1:end-2 end-2] );

% Growth derivative
dF = -( ( Ga.*Fa - Gb.*Fb )./ (ya - yb ) )';

dXdt = [dF; 0];

end % function