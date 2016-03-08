function status = solveroutput(t,~,flag)

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


% This prints an output during the course of solution

persistent START_TIME FINAL_TIME NUM_DOTS

if isequal(flag,'init')
    
    % Save the end time
    START_TIME = t(1);
    FINAL_TIME = t(end);
    
    % Print a header column
    fprintf('|%-10.8g %s %10.8g|\n|',START_TIME,repmat('.',1,28),FINAL_TIME)
    
    % No dots printed yet
    NUM_DOTS = 0;
    
elseif isempty(flag)
    
    % Print a number of dots depending on how far in time the solver has
    % reached
    
    % Calculate a percentage, rounded to the nearest 2 percent (50 dots
    % in total).
    progress = (max(t)-START_TIME)/(FINAL_TIME-START_TIME);
    
    % Calculate the total number of dots which should be printed up to
    % this time
    ndots = round(progress*50);
    
    if ndots > NUM_DOTS
        % Print the number of dots which are missing
        fprintf('%s',repmat('.',1,ndots-NUM_DOTS));
    end % if
    
    % Update the number of existing dots
    NUM_DOTS = ndots;
    
elseif isequal(flag,'done')
    
    % Check so that 50 dots were printed - otherwise, print the
    % remaining spaces
    if NUM_DOTS < 50
        fprintf('%s',repmat(' ',1,50-NUM_DOTS));
    end % if
    
    % Print the final |\n
    fprintf('|\n')
    
end % if elseif

% No cancelling
status = 0;

end % function