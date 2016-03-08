function O = solve(O)

% CATTube.gui_solve

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


% Disable 'Solve' button
set(O.gui.run.start,'Enable','off')

% Run the real solve method
solve@CAT(O);

% Activate the plot button and the overwrite checkbox
set(O.gui.run.plot,'Enable','on');
set(O.gui.run.plot_overwrite,'Enable','on');

end % function