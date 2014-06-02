function load(O,source)

% CAT.load
%
% Load a CAT object from a source (mat file or CAT/CATTube
% object) and import its settings into this instance of
% CAT/CATTube.
%
% EXAMPLE
%   CAT.load('saved-simulation.mat'); % Load from file
%
%   CAT.load(kitty); % Load from workspace CAT object
%
% If no source is given, the user is asked for a filename on
% the command line. If a source is given:
%  - a file is searched for any CAT* objects. The first one found is written into the current
%    instance of CAT.
%  - a CAT object is overwritten on top of the current object
%
% SEE ALSO
% CAT, CATTube, CAT.save

if nargin < 2 || isempty(source)
    % Allow an empty filename to get the filename from a GUI
    % window
    
    % Get the file name
    source = input('Enter the name of the .mat file to load: ','s');
    
end % if nargin

if ~isempty(source)
    
    if ischar(source) && exist(source,'file')
        
        % Get list of variables in the mat file - look for CAT/CATTube
        % objects
        
        S = whos('-file',source);
        
        Ssearch = regexp({S.class},'^CAT');
        
        if any([Ssearch{:}])
            % At least one CAT object found - load the first one.
            % This is maybe not ideal - but it works
            
            CATinput = S(find(Ssearch{:},1));
            Cname = CATinput.name;
            
            % Get this variable from the file
            CATinput = load(source,Cname);
            CATinput = CATinput.(Cname);
            
        else
            warning('CAT:load:nodata','The chosen file does not contain any useful data')
        end % if else
        
    elseif isa(source,'CAT')
        % A CAT object was given - use it!
        CATinput = source;
    else
        warning('CAT:load:badinput','The specified source is not valid')
    end % if
    
    % Overwrite the properties of the current object with
    % the values for those properties from the loaded
    % object
    fieldnames = properties(O);
    
    % TEMP FIX - turn warnings off to avoid warnings for empty fields
    warning('off','all')
    
    for i = 1:length(fieldnames)
        O.(fieldnames{i}) = CATinput.(fieldnames{i});
    end
    
    % TEMP FIX - turn warnings back on
    warning('on','all')
    
end % if - n warning needed here

end % function load