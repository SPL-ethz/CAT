function outputfilename = save(O,CATname)

if nargin < 2 || isempty(CATname)
    % CATname is the name of the object in the workspace
    CATname = inputname(1);
end % if

% Get the file name
[FileName,PathName] = uiputfile({'*.mat','MAT file (*.mat)';'*.m','Script (*.m)'},'Save .mat or .m file',pwd);
outputfilename = [PathName FileName];

if ~isempty(outputfilename) && ~isequal(outputfilename,[0 0])
    outputfilename = save@CAT(O,outputfilename,CATname);
    
    % Show a window
    helpdlg(sprintf('Saved file to %s\n',outputfilename),'CAT file saved');
    
    % Clear outputfilename if no output requested
    if nargout < 1
        clear outputfilename
    end % if
end % if

end % function