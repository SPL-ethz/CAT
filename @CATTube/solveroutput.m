function status = solveroutput(t,~,flag)

% This prints an output during the course of solution

persistent START_TIME FINAL_TIME HANDLE

if isequal(flag,'init')
    
    % Save the end time
    START_TIME = t(1);
    FINAL_TIME = t(end);
    
    % Create the wait bar
    HANDLE = waitbar(0,sprintf('Solving. Times: %g - %g',START_TIME,FINAL_TIME),...
    'CreateCancelBtn',@cancelintegration);
    
    % No cancelling
    status = 0;
    
elseif isempty(flag)
    
    % Print a number of dots depending on how far in time the solver has
    % reached
    
    % Check for the waitbar
    if exist('HANDLE','var') && ishandle(HANDLE)
        
        % Calculate a percentage, rounded to the nearest 2 percent (50 dots
        % in total).
        progress = (max(t)-START_TIME)/(FINAL_TIME-START_TIME);
        
        waitbar(progress,HANDLE)
        
        % Don't stop execution
        status = 0;
        
    elseif exist('HANDLE','var') && ~ishandle(HANDLE)
        % The handle has been closed/deleted by the cancel button - stop
        % execution
        status = 1;
    end % if
    
elseif isequal(flag,'done')
    
    % Close the wait bar - if it still exists
    if exist('HANDLE','var') && ishandle(HANDLE)
        delete(HANDLE)
    end % if
    
    % No cancelling
    status = 0;
    
end % if elseif

    function cancelintegration(varargin)

        delete(HANDLE)

    end % function

end % function