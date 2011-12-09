function [SolutionDists SolutionConc] = PBESolver ( ProblemDefinition )

% PBESOLVER
%
% Solve PBEs, like a boss

% Setup Solver - for constant y at the moment

X0 = [ProblemDefinition.init_dist.F ProblemDefinition.init_conc];

solvefun = str2func(ProblemDefinition.sol_method);

y = ProblemDefinition.init_dist.y;

solvefun = @(t,X) solvefun(t,X,y,ProblemDefinition.growthrate);

[t_out,X_out] = ode15s(solvefun , ProblemDefinition.sol_time , X0 );

% Create solution
for i = 1:length(t_out)
    SolutionDists(i) = Distribution( y , X_out(i,1:length(y)) );
end % for

SolutionConc = X_out(i,end);

end % function