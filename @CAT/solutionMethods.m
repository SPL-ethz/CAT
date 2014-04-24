function solmethods = solutionMethods(O)

% solutionMethods
%
% Print or return a list of the currently available solution methods

% Get all methods in CAT class
m = methods(O);

m = m( strncmp(m,'solver_',length('solver_')) );

for i = 1:length(m)
    solvname = m{i};
    m{i} = solvname(length('solver_')+1:end);
end % for

if nargout >= 1
    solmethods = m;
else
    disp(m)
end % if

end % function