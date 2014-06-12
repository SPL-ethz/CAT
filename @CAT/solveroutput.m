function status = solveroutput(t,~,flag)

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
    
    for i = 1:length(t)
        
        % Calculate a percentage, rounded to the nearest 2 percent (50 dots
        % in total).
        progress = (t(i)-START_TIME)/(FINAL_TIME-START_TIME);
        
        % Calculate the total number of dots which should be printed up to
        % this time
        ndots = round(progress*50);
        
        if ndots > NUM_DOTS
            % Print the number of dots which are missing
            fprintf('%s',repmat('.',1,ndots-NUM_DOTS));
        end % if
        
        % Update the number of existing dots
        NUM_DOTS = ndots;
        
    end % for
    
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