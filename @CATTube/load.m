function load(O,source)

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