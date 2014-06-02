function outputfilename = save(O,outputfilename,CATname)

% CAT.save
%
% This function saves the CAT instance in a mat file
%
% The name of the mat file can be given on the command line. If
% it is not given, a filename is created from the name of the
% object and the current date and time
%
% SEE ALSO
% CAT, CATTube, CAT.load

if nargin < 3 || isempty(CATname)
    CATname = inputname(1);
end % if

if nargin < 2 || isempty(outputfilename)
    if isempty(CATname)
        % CATname can still be empty if no object name could be
        % found (result of calculation on command line)
        CATname = 'CAT';
    end % if
    outputfilename = strcat(CATname,'_',datestr(now,'YYYYmmdd_HHMMSS'),'.mat');
end % if

% Check the path and add .mat if necessary
if ~regexp(outputfilename,'\.mat$')
    outputfilename = [outputfilename '.mat'];
end % if

superCAT.(CATname) = O; %#ok<STRNU>
save(outputfilename,'-struct','superCAT')

% Clear outputfilename if no output requested
if nargout < 1
    clear outputfilename
    % Acknowledge saving
    fprintf('Saved file to %s\n',outputfilename);
end % if

end % function