function setDefaults(O)

% setDefaults
%
% set sensible defaults for CAT method

O.init_dist = Distribution(linspace(1,500,200),{'normal',100,20});
O.init_conc = 1; % saturated
O.solubility = @(T,xm) 1;
O.Tprofile = @(t) 25*ones(size(t));
O.ASprofile = @(t) 0*t;
O.growthrate = @(S,T,y) ones(size(y));
O.nucleationrate = @(S,T,m) 0;
O.init_seed = 1;
O.init_massmedium = 1000;
O.rhoc = 1e-12;
O.kv = 1;
O.sol_method = 'centraldifference';
O.sol_time = [0 100];

end