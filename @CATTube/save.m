function outputfilename = save(O,CATname)

if nargin < 2 || isempty(CATname)
    % CATname is the name of the object in the workspace
    CATname = inputname(1);
end % if


% Get the file name
[FileName,PathName] = uiputfile('*.mat','Save .mat file',pwd);
outputfilename = [PathName FileName];

if ~isempty(outputfilename) && ~isequal(outputfilename,[0 0])
    outputfilename = save@CAT(O,outputfilename,CATname);
end % if

% Show a window
helpdlg(sprintf('Saved file to %s\n',outputfilename),'CAT File saved');

% Clear outputfilename if no output requested
if nargout < 1
    clear outputfilename
end % if

end % function