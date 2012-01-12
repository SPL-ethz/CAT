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
        
        %% Results
        
        % Vector of actual times returned by solver
        calc_time
        
        % Distributions for each time step
        calc_dist
        
        % Concentrations over time
        calc_conc
        
        % Temperature over time
        calc_temp
        
        % Volume (of reactor content) over time
        calc_volume
        
        % Method to use - default to central difference
        sol_method = 'centraldifference'
        
        sol_options = {};
        
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
        coolingrate = @(t) 0;
        
        % anti-solvent addition rate
        % Defines the anti-solvent addition rate dV/dt.         
        ASadditionrate = 0;
        
    end % properties
    
    methods
        
        %% Method ProblemDefinition (constructor)
        
        function O = ProblemDefinition(init_dist,init_conc,sol_time,...
                sol_method,growthrate, nucleationrate, coolingrate, ASadditionrate)
            
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
            
            if nargin > 6 && ~isempty(coolingrate)
                O.coolingrate = coolingrate;
            end % if 
            
            if nargin > 7 && ~isempty(ASadditionrate)
                O.ASadditionrate = ASadditionrate;
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
        
        function set.sol_options(O,value)
            
            % SET.SOL_OPTIONS
            %
            % Check the defined OPTIONS. Should be a cell array
            %
            % Probable more checks should be carried out at this point in
            % the future
            
            if iscell(value)
                O.sol_options = value;
            else
                warning('ProblemDefinition:SetSol_Options:WrongType',...
                    'The property sol_options should be a cell');
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
        
        %% Method set.coolingrate
        
        function set.coolingrate(O,value)
            
            % SET.GROWTHRATE
            %
            % Check the cooling rate: either it is a scalar (for cooling it
            % should be less than zero) or a function handle
            % arguments: t (scalar)
            % output: dTdt (scalar)
            
            if strcmp(class(value),'function_handle')                
                % Check the number of inputs
                if nargin(value) ~=1                    
                    warning('Distribution:setcoolingrate:Wrongnargin',...
                            'The cooling rate must have a single input parameter (time).');
                else
                    % The growth rate function is OK, set it
                    O.coolingrate = value;
                end % if else 
            elseif isscalar(value)
                O.coolingrate = @(t) double(value);                
            else
                warning('Distribution:setcoolingrate:Wrongtype',...
                    'The cooling rate must be defined as a function_handle with 1 input or a scalar for a constant cooling rate.');
            end %if
            
        end % function
        
        %% Method massbal
        
        function PDma = massbal(O)

            PDma = abs(moments(O.calc_dist,3,1)*O.rhoc*O.kv+O.init_conc-...
                (moments(O.calc_dist,3)*O.rhoc*O.kv+O.calc_conc'))./...
                (moments(O.calc_dist,3,1)*O.rhoc*O.kv+O.init_conc);
            
        end % function
        
                %% Method plot
        
        function PDpl = plot(O,plotwhat,varargin)
            
            % Plot function for the results
            %
            % Use plot(PD,plotwhat) to plot the results of the
            % simulation. plotwhat is a string that defines what exactly
            % should be plotted. Possible input:
            % 'results'         -   plot everything
            % 'distributions'   -   plot distributions
            % 'distoverlap'     -   only plot overlapping distributions (2D)
            % 'dist3D'          -   only plot 3D surf plot of distributions
            % 'cumprop'         -   plot cumulative properties (moments)
            % 'process'         -   plot process variables (T, conc)
            % 
            % and any combination thereof.
            %
            % PLOT returns the handles to the plot objects created. (If
            % several plots are created the sequence of handles is the
            % following: PSDs overlapping, PSDs3D, cumulative properties
            % (moments), process variables (temperature, concentration)
            %
            % Graphs can currently not be plotted in existing figures !!
            %
                
                PDpl = [];
                
            if (~isempty(find(strcmp(plotwhat,'distributions'))) ...
                    || ~isempty(find(strcmp(plotwhat,'distoverlap'))) ...
                    || ~isempty(find(strcmp(plotwhat,'results'))))
                %% currently not active
%                 % Check if Parent axes are already defined
%                 useaxpos = find(strcmp(varargin,'Parent'));
% 
%                 if ~isempty(useaxpos) && ishandle(varargin{useaxpos+1})
%                     % Define this axes as the one to use
%                     Fax = varargin{useaxpos+1};
%                     % Remove this parent command from varargin, is added again
%                     % later
%                     varargin(useaxpos+(0:1)) = [];
%                 else
%%
                    FFig = figure(11);
                    set(gcf,'numbertitle','off','name','PSDs (overlapping)')
                    Fax(1) = subplot(1,2,1);
                    Fax(2) = subplot(1,2,2);
                    xlabel(Fax(1),'Mean Char. Length')
                    xlabel(Fax(2),'Mean Char. Length')
                    ylabel(Fax(1),'Number Distribution')
                    ylabel(Fax(2),'Normalized Volume Distribution')
%                 end % if

                % Handles for plots
                PDpl_local = zeros(1,length(O.calc_dist));

                hold(Fax(1),'all')
                hold(Fax(2),'all')

                % Plot every distribution
                for i = 1:length(O.calc_dist)
                    PDpl_local(i) = plot(...
                        O.calc_dist(i).y,O.calc_dist(i).F,...
                        'Parent',Fax(1),'DisplayName',['Dist ' num2str(i)],...
                        varargin{:});
                    
                    PDpl_local(length(O.calc_dist)+i) = plot(...
                        O.calc_dist(i).y,O.calc_dist(i).F.*O.calc_dist(i).y.^3./moments(O.calc_dist(i),3),...
                        'Parent',Fax(2),'DisplayName',['Dist ' num2str(i)],...
                        varargin{:});
                end % for
                
                PDpl = [PDpl PDpl_local];
            end % if
            
  
            % 3D plot of distributions over time
            if (~isempty(find(strcmp(plotwhat,'distributions'))) ...
                    || ~isempty(find(strcmp(plotwhat,'dist3D'))) ...
                    || ~isempty(find(strcmp(plotwhat,'results')))...
                    && ~strcmp(O.sol_method,'movingpivot'))
                
                for i = 1:length(O.calc_dist)
                    Fmat(:,i) = O.calc_dist(i).F;
                end % for
                
                figure(12)
                set(gcf,'numbertitle','off','name','PSDs (3D time evolution)')
                
                % Handles for plots
                PDpl_local = zeros(1,2);

                subplot(1,2,1)
                PDpl_local(1) = surf(O.calc_time(:),O.calc_dist(1).y(:),Fmat,varargin{:});
                ylabel('Mean Char. Length')
                xlabel('Time')
                zlabel('Number Distribution')

                subplot(1,2,2)
                PDpl_local(2) = surf(O.calc_time,O.calc_dist(1).y,...
                    Fmat.*repmat(O.calc_dist(1).y(:).^3,1,size(O.calc_time))...
                    ./repmat(moments(O.calc_dist,3),length(O.calc_dist(1).y),1),...
                    varargin{:});
                ylabel('Mean Char. Length')
                xlabel('Time')
                zlabel('Normalized Volume Distribution')

                PDpl = [PDpl PDpl_local];
            end % if
            
            % Cumulative Properties
            if (~isempty(find(strcmp(plotwhat,'results'))) || ...
                ~isempty(find(strcmp(plotwhat,'cumprop'))))
                
                figure(21)
                set(gcf,'numbertitle','off','name','PSD cumulative properties')  

                % Handles for plots
                PDpl_local = zeros(1,3);
                
                subplot(3,1,1)
                PDpl_local(1) = plot(O.calc_time,moments(O.calc_dist,0));
                title(strcat('Setup: J = ',num2str(exist('nucleationrate'))))
                ylabel('0^{th} moment')
                if ~exist('nucleationrate')
                    ylim([(1-0.015)*moments(O.calc_dist(1),0) 1.015*moments(O.calc_dist(1),0)])
                    set(gca,'ytick',[(1-0.01)*moments(O.calc_dist(1),0) moments(O.calc_dist(1),0) 1.01*moments(O.calc_dist(1),0)],'yticklabel',{'-1%' '0%' '+1%'});
                end % if

                subplot(3,1,2)
                PDpl_local(2) = plot(O.calc_time,moments(O.calc_dist,3));
                ylabel('3^{rd} moment')
                
                subplot(3,1,3)
                PDpl_local(3) = plot(O.calc_time,moments(O.calc_dist,4)./moments(O.calc_dist,3));
                ylabel('Weight average length')
                xlabel('Time')
                
                PDpl = [PDpl PDpl_local];
            end % if
            
            % Process Variables
            if (~isempty(find(strcmp(plotwhat,'results'))) || ...
                ~isempty(find(strcmp(plotwhat,'process'))))
            
                figure(31)
                set(gcf,'numbertitle','off','name','Process Variables')
            
            % Handles for plots
                PDpl_local = zeros(1,1);
            
                if ~exist('O.calc_cstar')
                    
                    subplot(1,3,1);                    
                    PDpl_local = plot(O.calc_time,O.calc_conc);
                    xlabel('Time')
                    ylabel('Concentration')
                    PDpl = [PDpl PDpl_local];
                    
                    subplot(1,3,2);                    
                    PDpl_local = plot(O.calc_time,O.calc_temp);
                    xlabel('Time')
                    ylabel('Temperature')
                    PDpl = [PDpl PDpl_local];
                    
                    subplot(1,3,3);                    
                    PDpl_local = plot(O.calc_time,O.calc_volume);
                    xlabel('Time')
                    ylabel('Volume of reactor content')
                    PDpl = [PDpl PDpl_local];
                else
                    warning('ProblemDefinition:PlotCumProp:Inexistent',...
                    'Concentration profile output is missing');
                end % if
            end % if
            
        end % function
        
    end % methods
    
end % classdef