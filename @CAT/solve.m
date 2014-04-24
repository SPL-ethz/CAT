function solve(O)
%% CAT.solve
% solve is a method of the CAT class. It integrates the problem according to
% the settings of the CAT object and fills in the calc-properties of the
% same object. The method does not accept any additional inputs, all options 
% should therefore be set via the property fields of the object.
% Solve is structured into two nested parts: An upper layer that handles
% discontinuities in the supplied AS or T profiles and a lower layer (the
% solvers).
%% Solve
O.calc_dist = Distribution;
if ~isempty(O.init_seed)
    O.init_dist.mass = [O.init_seed O.kv O.rhoc O.init_massmedium];
end

if ~isempty(O.tNodes)
    
    % Save for later
    sol_time = O.sol_time;
    
    for i = 2:length(O.tNodes) % make sure you hit the different nodes of the non-smooth profiles
        
        % Cut out the piece we want to look at currently
        O.sol_time = [O.tNodes(i-1) sol_time(sol_time>O.tNodes(i-1) & sol_time<O.tNodes(i)) O.tNodes(i)];
        
        % Solve for the current piece
        try
            % Simply run the method corresponding to the chosen solution method
            % prefixed with solver_
            O.(['solver_' O.sol_method]);
        catch ME
            error('solve:tryconsttemp:PBESolverfail',...
            'Solver failed to integrate your problem. Message: %s',ME.message)
        end
        
        % Save the solution for the current piece
        calc_time(end+1:end+length(O.calc_time)) = O.calc_time;
        calc_dist(end+1:end+length(O.calc_dist)) = O.calc_dist;
        calc_conc(end+1:end+length(O.calc_conc)) = O.calc_conc; 
        
        % Set new initial distribution and initial concentration
        O.init_dist = O.calc_dist(end);
        O.init_conc = O.calc_conc(end);

    end % for

    % Put temporary solution into the right place
    O.calc_time = calc_time;
    O.calc_dist = calc_dist;
    O.calc_conc = calc_conc;
    % Reset times and initial distribution, concentration
    O.sol_time = sol_time;
    O.init_dist = calc_dist(1);
    O.init_calc = calc_conc(1);
    
else    
    
    try
        
        % Simply run the method corresponding to the chosen solution method
        % prefixed with solver_
        O.(['solver_' O.sol_method]);
        
    catch ME
        error('solve:tryconsttemp:PBESolverfail',...
            'Solver failed to integrate your problem. Message: %s',ME.message)
    end
end

O.init_conc = O.calc_conc(1);
O.init_dist = O.calc_dist(1);

if length(O.sol_time)>2
    [~,I] = intersect(O.calc_time,O.sol_time);
    O.calc_time = O.calc_time(I);
    O.calc_dist = O.calc_dist(I);
    O.calc_conc = O.calc_conc(I);
end

%% Check Mass balance
if any(O.massbal > 5)
   warning('ProfileManager:massbalcheck:largeerror',...
                    'Your mass balance error is unusually large (%4.2f%%). Check validity of equations and consider increasing the number of bins.',max(O.massbal)); 
end

end