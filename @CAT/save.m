function outputfilename = save(O,outputfilename,CATname)

% CAT.save
%
% This function saves the CAT instance in a mat or m file. The default is
% to save the object to a mat file.
%
% If a filename ending with .m is specified, a script is created which
% allows the user to recreate the current object.
%
% The name of the output file can be given on the command line. If
% it is not given, a filename is created from the name of the
% object and the current date and time.
%
% SEE ALSO
% CAT, CATTube, CAT.load, CAT.saveSource

if nargin < 3 || isempty(CATname)
    CATname = inputname(1);
end % if

if nargin < 2 || isempty(outputfilename)
    if isempty(CATname)
        % CATname can still be empty if no object name could be
        % found (result of calculation on command line)
        CATname = 'kitty';
    end % if
    outputfilename = strcat(CATname,'_',datestr(now,'YYYYmmdd_HHMMSS'),'.mat');
end % if

% Check for .m file ending - this means that user wants to write a script
% file
if regexp(outputfilename,'\.m$')
    saveSource(O,outputfilename,CATname);
else
    % Default is to assume a .mat-file is to be written
    
    % Check the path and add .mat if necessary
    if ~regexp(outputfilename,'\.mat$')
        outputfilename = [outputfilename '.mat'];
    end % if
    
    % Save the mat file
    saveMAT(O,outputfilename,CATname);
    
end % if

% Clear outputfilename if no output requested
if nargout < 1
    % Acknowledge saving
    fprintf('Saved file to %s\n',outputfilename);
    clear outputfilename
end % if

end % function

%% Function saveMAT

function saveMAT(O,outputfilename,CATname)

% saveMAT
%
% Auxiliary function which takes care of saving the CAT object to a .mat
% file

superCAT.(CATname) = O; %#ok<STRNU>
save(outputfilename,'-struct','superCAT')

end % saveMAT

%% Function saveSource

function saveSource(O,outputfilename,CATname)

% saveSource
%
% Auxiliary function which takes care of saving the CAT object to a script
% file (writes source)

% Open the output file for writing - discard 
fid = fopen(outputfilename,'w+');

if fid < 3 % fid needs to be at least 3 if the file was opened successfully
    error('CAT:save:failedopenoutput','The chosen file could not be opened for writing. Do you have the correct permissions?');
end % if

% Output the header
fprintf(fid,[...
    '%%%% CAT settings\n' ...
    '%% Crystallization Analysis Toolbox\n' ...
    '%% File created automatically on %s\n' ...
    '%% \n'...
    '%% CAT Toolbox: http://www.ipe.ethz.ch/laboratories/spl/cat\n\n']...
    ,datestr(now));

% Initialise the CAT object
fprintf(fid,'%s = CAT;\n\n',CATname);

% Get list of properties, modify it to exlude things which are not be
% written to a file
fieldnames = properties(O);
fieldnames(strcmp(fieldnames,'gui')) = [];
fieldnames(~cellfun(@isempty,strfind(fieldnames,'calc'))) = [];

% To make the output:
%   1. Check for multidimensional CAT objects - if this is the case, run
%   the output section for properties several times and update CATname in
%   between each run
%   2. Make a new cell for the next entry
%
% For each property
%   1. Print the property help (this includes the property name) as a
%   comment
%   2. Add the output property

for iv = 1:length(O)
    
    % Print cell header
    if length(O) > 1
        fprintf(fid,'%%%% Properties (%i of %i)\n\n',iv,length(O));
        CATnamevec = sprintf('%s(%i)',CATname,iv);
    else
        fprintf(fid,'%%%% Properties\n\n');
        CATnamevec = CATname;
    end % if
    
    for i = 1:length(fieldnames)
        
        % Get help string
        helpstr = help(strcat('CAT.',fieldnames{i}));
        
        % Remove the last line break and add comment signs at the beginning of
        % each line
        helpstr = strrep(['%' helpstr(1:end-1)],char(10),[10 '%']);
        
        % Print the output for this property
        fprintf(fid,'%s\n%s.%s = %s;\n\n',...
            helpstr,CATnamevec,fieldnames{i},data2str(O(iv).(fieldnames{i}))...
            );
        
    end % for
    
end % for

% Print solve footer
fprintf(fid,'%%%% Solve\n\n%s.solve;\n',CATname);

% Close the file
fclose(fid);

end % function saveSource