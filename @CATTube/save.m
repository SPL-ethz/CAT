function outputfilename = save(O,outputfilename,CATname)

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


if nargin < 3 || isempty(CATname)
    % CATname is the name of the object in the workspace
    CATname = inputname(1);
end % if

if nargin < 2 || isempty(outputfilename)
    % Get the file name
    [FileName,PathName] = uiputfile({'*.mat','MAT file (*.mat)';'*.m','Script (*.m)'},'Save .mat or .m file',pwd);
    outputfilename = [PathName FileName];
end % if

if ~isempty(outputfilename) && ~isequal(outputfilename,[0 0])
    
    % Want to save data as CAT class, not as CATTube. Some info will get
    % lost this way but probably wise to have only one class for storage.
    OCat = CAT;
    OCat.load(O); % reads content of current CATTube object into new CAT instance
    
    outputfilename = OCat.save(outputfilename,CATname,O.gui.source);
    
    % Show a window
    helpdlg(sprintf('Saved file to %s\n',outputfilename),'CAT file saved');
    
    % Clear outputfilename if no output requested
    if nargout < 1
        clear outputfilename
    end % if
end % if

end % function