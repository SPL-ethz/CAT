classdef ProblemDefinition < handle
    
    % ProblemDefinition
    
    properties
        
        % Initial distribution (Distribution object, defines size and
        % distribution
        init_dist = Distribution
        
        % Initial concentration
        init_conc = 1;
        
        % Initial temperature
        init_temp = 298;
        
        % Initial volume
        init_volume = 1;
        
        % Solution time vector
        sol_time
        
        % Crystal density
        rhoc = 1000
        
        % Shape factor
        kv = 1
        
        % Seed Mass
        seed_mass = 1; 
        
        % Method to use - default to central difference
        sol_method = 'centraldifference'
        
        % Growth rate function
        % This function should be called as growthrate(c,y) where c is the current
        % concentration and y is the size. It should return a vector the
        % same size as y
        growthrate = @(c,T,y) ones(size(y))
        
        % Nucleation rate function
        % This function should be called as nucleationrate(c,T) where c is the current
        % concentration and T is the temperature. It returns a scalar.
        nucleationrate = @(c,T) 0
        
        % cooling rate
        % Defines the cooling rate dT/dt. 
        
        coolingrate = 0;
        
    end % properties
    
    methods
        
        %% Method ProblemDefinition (constructor)
        
        function O = ProblemDefinition(init_dist,init_conc,sol_time,...
                sol_method,growthrate, nucleationrate)
            
            % PROBLEMDEFINITION
            %
            % Constructor function for ProblemDefinition class. Called as:
            %	PD = ProblemDefinition(init_dist,init_conc,sol_time,...
            %      sol_method,growthrate)
            %
            % Leave any of the variables empty to not assign them.
            %
            % The variables which are given are assigned directly to the
            % corresponding property. All of them can also be assigned
            % later.
            
            % Set values, if they are given
            
            if nargin > 0 && ~isempty(init_dist)
                O.init_dist = init_dist;
            end % if
            
            if nargin > 1 && ~isempty(init_conc)
                O.init_conc = init_conc;
            end % if
            
            if nargin > 2 && ~isempty(sol_time)
                O.sol_time = sol_time;
            end % if
            
            if nargin > 3 && ~isempty(sol_method)
                O.sol_method = sol_method;
            end % if
            
            if nargin > 4 && ~isempty(growthrate)
                O.growthrate = growthrate;
            end % if
            
            if nargin > 5 && ~isempty(nucleationrate)
                O.nucleationrate = nucleationrate;
            end % if    
        end % function
        
        %% Method set.init_dist
        
        function set.init_dist(O,value)
            
            % SET.INIT_DIST
            %
            % Check the initial distribution, it must be a distribution
            % class object
            
            if ~strcmp(class(value),'Distribution')
                warning('ProblemDefinition:SetInit_Dist:WrongType',...
                    'The init_dist property must be a Distribution object');
            else
                O.init_dist = value;
            end % if else
            
        end % function
        
        %% Method set.init_conc
        
        function set.init_conc(O,value)
            
            % SET.INIT_CONC
            %
            % Check the initial concentration, it must be a positive,
            % finite scalar (can be zero)
            
            if isscalar(value) && value >= 0 && isfinite(value)
                O.init_conc = value;
            else
                warning('ProblemDefinition:SetInit_Conc:WrongType',...
                    'The init_conc property must be a positive, finite scalar (may be zero)');
            end % if else
            
        end % function
        
        %% Method set.sol_time
        
        function set.sol_time(O,value)
            
            % SET.SOL_TIME
            %
            % Check the solution time vector. This should be a vector of
            % monotonically increasing values
            
            if length(value) > 1 && isvector(value) && ~any(diff(value)<=0)
                O.sol_time = value;
            else
                warning('ProblemDefinition:SetSol_Time:WrongValue',...
                    'The property sol_time must be a vector of monotonically increasing values');
            end % if else
            
        end % function
        
        %% Method set.method
        
        function set.sol_method(O,value)
            
            % SET.SOL_METHOD
            %
            % Check the defined solution method. Should be a string
            %
            % Probable more checks should be carried out at this point in
            % the future
            
            if ischar(value)
                O.sol_method = value;
            else
                warning('ProblemDefinition:SetSol_Method:WrongType',...
                    'The property sol_method should be a string');
            end % if else
            
        end % function
                
        %% Method set.growthrate
        
        function set.growthrate(O,value)
            
            % SET.GROWTHRATE
            %
            % Check the growth rate: should be a function handle, accept 3
            % arguments: c (scalar), T (temperature), and y (vector). The output should be
            % the same size as y
            
            if strcmp(class(value),'function_handle')
                
                % Check the number of inputs
                if nargin(value) == 3
                    
                    % Check the output using 2 example values
                    out = value(1.1,1,linspace(0.1,1,10));
                    
                    % Check size
                    if any( size(out) ~= [1 10] )
                        % Size of output wrong
                        warning('Distribution:setgrowthrate:Wrongsize',...
                            'The growth rate function returns a vector which is not the same size as the input vector');
                    else
                        % The growth rate function is OK, set it
                        O.growthrate = value;
                    end % if else
                    
                else
                    warning('Distribution:setgrowthrate:Wrongnargin',...
                        'The growth rate function must have 3 input arguments (concentration, temperature, sizes)');
                end % if
                
            else
                warning('Distribution:setgrowthrate:Wrongtype',...
                    'The growth rate must be defined as a function');
            end %if
            
        end % function
        
    end % methods
    
end % classdef