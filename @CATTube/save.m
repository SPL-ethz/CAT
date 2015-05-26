function outputfilename = save(O,outputfilename,CATname)

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
    O.load(O); % reads content of current CATTube object into new CAT instance
    
    outputfilename = OCat.save(outputfilename,CATname,O.gui.source);
    
    % Show a window
    helpdlg(sprintf('Saved file to %s\n',outputfilename),'CAT file saved');
    
    % Clear outputfilename if no output requested
    if nargout < 1
        clear outputfilename
    end % if
end % if

end % function