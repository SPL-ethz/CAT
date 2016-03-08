function load(O,source)

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


if nargin < 2 || isempty(source)
    % Allow an empty source to get the filename from a GUI
    % window
    
    % Get the file name
    [FileName,PathName] = uigetfile('*.mat','Load .mat file',pwd);
    
    source = [PathName FileName];
    
end % if nargin

if ~isequal(source,[0 0])
    % When 'cancel' button is pressed, output of GUI function
    % is [0 0] - don't send this to load function, only if file
    % chosen
    
    % Send to CAT version of load
    load@CAT(O,source)
end % if

end % function