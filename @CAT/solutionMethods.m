function solmethods = solutionMethods(O)

% solutionMethods
%
% Print or return a list of the currently available solution methods

% Copyright 2015-2016 David Ochsenbein
% Copyright 2012-2014 David Ochsenbein, Martin Iggland
% 
% This file is part of CAT.
% 
% CAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% 
% CAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


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