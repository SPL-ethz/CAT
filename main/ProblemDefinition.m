classdef ProblemDefinition < handle
    
    % ProblemDefinition
    
    properties
        
        % Initial distribution (Distribution object, defines size and
        % distribution
        init_dist = Distribution
        
        % Initial concentration
        init_conc = 1;
        
        % Solubility
        solubility = @(T,xm) 0.5;
        
        % Temperature profile
        Tprofile = @(t) 298*ones(size(t));
        
        % (Anti)solvent added mass profile
        ASprofile = @(t) 0*t;

        % nodes for non-smooth input profiles
        tNodes = [];
        
        % Seed mass
        init_seed = 1;
        
        % Initial mass of solvent + antisolvent
        init_massmedium = 1;
        
        % Solution time vector
        sol_time
        
        % Crystal density [g/micron^3]
        rhoc = 1e-12
        
        % Shape factor
        kv = 1
        
        %% Results
        
        % Vector of actual times returned by solver
        calc_time
        
        % Distributions for each time step
        calc_dist
        
        % Concentrations over time [g solute / g total solvents]
        calc_conc
        
        % Method to use - default to central difference
        sol_method = 'centraldifference'
        
        % Solver options -  default none
        sol_options = {};
        
        % Growth rate function This function should be called as
        % growthrate(S,T,y) where S is the current supersaturation, T is the
        % temperature and y is the size. It should return a vector the same
        % size as y
        growthrate = @(S,T,y) ones(size(y))
        
        % Nucleation rate function
        % This function should be called as nucleationrate(c,T,m) where S is the current
        % supersaturation, T is the temperature. Optionally, the user can
        % specificy that the nucleation rate depends on a moment m of the
        % passed distribution F
        nucleationrate = @(S,T,m,F) 0
        
  
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
        
        %% Method set.init_seed
        
        function set.init_seed(O,value)
            
            % SET.INIT_SEED
            %
            % Set mass of seeds
            
            if isscalar(value) && ~isnan(value) && value>0 && ~isinf(value)
                O.init_seed = value;
            else
                warning('ProblemDefinition:SetInit_Seed:WrongType',...
                    'The init_seed property must be a positive scalar');
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
                
            elseif ischar(value) && strcmpi(value,'sat') % solution is saturated in the beginning
                O.init_conc = O.solubility(O.Tprofile(0));
            else
                warning('ProblemDefinition:SetInit_Conc:WrongType',...
                    'The init_conc property must be a positive, finite scalar (may be zero) or the string ''sat''');
            end % if else
            
        end % function
        
        %% Method set.solubility
        
        function set.solubility(O,value)
            
            % SET.solubility
            %
            %  Setter method for solubility Must be a function handle with
            %  1 or 2 inputs
            if isa(value,'function_handle') && nargin(value)==1

                O.solubility = @(T,xm) value(T); % we assume the only input argument is Temperature
                
            elseif isa(value,'function_handle') && nargin(value)==2

                O.solubility = value;


            else
                warning('ProblemDefinition:SetASprofile:WrongType',...
                    'The ASprofile property must be a positive, finite matrix (may be zero) or a function handle with one input');
                
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
        
        %% Method set.sol_method
        
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
        
        %% Method set.sol_options
        
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
        
        function set.Tprofile(O,value)
            
            % SET.Tprofile
            %
            % Check the Temperature profile. It must be a matrix with
            % positive, finite elements. The first row indicates the times
            % of the nodes whereas the second row indicates Temp's

            if ismatrix(value) && length(value(:,1))==2 && all(isfinite(value(:)))
               
                if value(1,end)<O.sol_time(end)
                    value = [value [0;0]];
                    value(1,end) = O.sol_time(end);
                    value(2,end) = value(2,end-1);
                end
                O.Tprofile = @(t) interp1(value(1,:),value(2,:),t); %
                O.tNodes = unique([O.tNodes value(1,:)]);
                
                


            elseif isa(value,'function_handle') && nargin(value)==1
                O.Tprofile = value;

            else
                warning('ProblemDefinition:SetTprofile:WrongType',...
                    'The Tprofile property must be a positive, finite matrix (may be zero) or a function handle with one input');
                
            end % if else
            
            
        end % function
        
        
        function set.ASprofile(O,value)
            
            % SET.ASprofile
            %
            % Check the time profile (for added AS profiles). It must be
            % a strictly increasing(!), positive vector
            if (ismatrix(value) && length(value(:,1))==2 && all(value >= 0) && all(isfinite(value)) && all(diff(value(2,:))>0)) || (isa(value,'function_handle') && nargin(value)==1)

                O.ASprofile = @(t) interp1(value(1,:),value(2,:),t);

                O.tNodes = unique([O.tNodes value(1,:)]);

            else
                warning('ProblemDefinition:SetASprofile:WrongType',...
                    'The ASprofile property must be a positive, finite matrix (may be zero) or a function handle with one input');
                
            end % if else
            
        end % function

        function set.tNodes(O,value)
            
            % SET.tNodes
            %
            % Set time nodes (make sure integrator covers them correctly)
            if isvector(value) && all(value>=0)

                O.tNodes = value;

            else
                warning('ProblemDefinition:SettNodes:WrongType',...
                    'The tNodes property must be a non-negative vector');
                
            end % if else
            
        end % function
        
        %% Method set.growthrate
        
        function set.growthrate(O,value)
            
            % SET.GROWTHRATE
            %
            % Check the growth rate: should be a function handle, accept 3
            % arguments: S (scalar), T (temperature), and y (vector). The output should be
            % the same size as y
            
            if isa(value,'function_handle')
                
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
                        'The growth rate function must have 3 input arguments (supersaturation, temperature, sizes)');
                end % if
                
            else % not a function handle
                warning('Distribution:setgrowthrate:Wrongtype',...
                    'The growth rate must be defined as a function');
            end %if
            
        end % function

        %% Method set.growthrate
        
        function set.nucleationrate(O,value)
            
            % SET.nucleationrate
            %
            % Check the nucleationrate rate: should be a function handle, accept 4
            % arguments: S (supersaturation), T (temperature), and m (moment), F (distribution). The output should be
            % a scalar
            
            if isa(value,'function_handle')
                
                % Check the number of inputs
                if nargin(value) == 4
                    
                    O.nucleationrate = value;
                    
                elseif nargin(value) == 2
                    
                    O.nucleationrate = @(S,T,~,~) value(S,T);
                end
                
            else % not a function handle
                warning('Distribution:setgrowthrate:Wrongtype',...
                    'The growth rate must be defined as a function');
            end %if
            
        end % function


        %% Method massbal
        
        function PDma = massbal(O)
            mscalc = O.init_massmedium+O.ASprofile(O.calc_time)-O.ASprofile(0);
            mass_solute = O.calc_conc(:).*mscalc(:);   % total mass of solute
            m3 = moments(O.calc_dist,3);                            % third moment over time
            mass_crystals = O.rhoc*O.kv*mscalc(:).*m3(:);    % total mass of crystals
            
            PDma = 100*((mass_solute + mass_crystals)/(mass_solute(1)+mass_crystals(1))-1);     % mass balance error   
        end % function
        
                %% Method plot
        
        function PDpl = plot(O,plotwhat,varargin)
            
            % Plot function for the results
            %
            % Use plot(PD,plotwhat) to plot the results of the
            % simulation. plotwhat is a string that defines what exactly
            % should be plotted. Possible input:
            % 'results'         -   plot everything
            % 'detailed_results'-   plot everything and more
            %
            % plotted using results
            % 'distributions'   -   plot distributions
            % 'distoverlap'     -   only plot overlapping distributions (2D)
            % 'dist3D'          -   only plot 3D surf plot of distributions
            % 'cumprop'         -   plot cumulative properties (moments)
            % 'process'         -   plot process variables (T, conc)
            %
            % additionally plotted in detailed mode
            % 'moments'         -   plots of the first four moments
            % 'integration'     -   details from the integration
            % (massbalance over time,...)
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
                Tcalc = O.Tprofile(O.calc_time);
                mscalc = O.init_massmedium+O.ASprofile(O.calc_time)-O.ASprofile(0);
                
            if (~isempty(find(strcmp(plotwhat,'distributions'))) ...
                    || ~isempty(find(strcmp(plotwhat,'distoverlap'))) ...
                    || ~isempty(find(strcmp(plotwhat,'detailed_results')))...
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
            % Note this feature is diasbled for moving pivot
            if (~isempty(find(strcmp(plotwhat,'distributions'))) ...
                    || ~isempty(find(strcmp(plotwhat,'dist3D'))) ...
                    || ~isempty(find(strcmp(plotwhat,'results')))...
                    || ~isempty(find(strcmp(plotwhat,'detailed_results')))...
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
                    Fmat.*repmat(O.calc_dist(1).y(:).^3,1,length(O.calc_time))...
                    ./repmat(moments(O.calc_dist,3),length(O.calc_dist(1).y),1),...
                    varargin{:});
                ylabel('Mean Char. Length')
                xlabel('Time')
                zlabel('Normalized Volume Distribution')

                PDpl = [PDpl PDpl_local];
                
            elseif (~isempty(find(strcmp(plotwhat,'distributions'))) ...
                    || ~isempty(find(strcmp(plotwhat,'dist3D'))) ...
                    || ~isempty(find(strcmp(plotwhat,'results')))...
                    || ~isempty(find(strcmp(plotwhat,'detailed_results')))...
                    && strcmp(O.sol_method,'movingpivot'))
                
                figure(12)
                set(gcf,'numbertitle','off','name','PSDs (3D time evolution)')
                
                % Handles for plots
                PDpl_local = zeros(1,2);

                subplot(1,2,1)
                for i = 1:length(O.calc_time)
                    PDpl_local(i) = plot3(repmat(O.calc_time(i),size(O.calc_dist(i).y)),O.calc_dist(i).y(:),O.calc_dist(i).F,varargin{:});
                    hold on
                end
                grid on
                hold off
                ylabel('Mean Char. Length')
                xlabel('Time')
                zlabel('Number Distribution')

                subplot(1,2,2)
                for i = 1:length(O.calc_time)
                    PDpl_local(i+length(O.calc_time)) = plot3(repmat(O.calc_time(i),size(O.calc_dist(i).y)),O.calc_dist(i).y,...
                    O.calc_dist(i).F(:).*O.calc_dist(i).y(:).^3./...
                    moments(O.calc_dist,3,i),...
                    varargin{:});
                hold on
                end
                grid on
                hold off
                ylabel('Mean Char. Length')
                xlabel('Time')
                zlabel('Normalized Volume Distribution')
                set(PDpl_local,'linewidth',1.5,'color','k')
                PDpl = [PDpl PDpl_local];
                
            end % if
            
            % Cumulative Properties
            if (~isempty(find(strcmp(plotwhat,'results'))) || ...
                ~isempty(find(strcmp(plotwhat,'detailed_results'))) || ...
                ~isempty(find(strcmp(plotwhat,'cumprop'))))
                
                figure(21)
                set(gcf,'numbertitle','off','name','PSD cumulative properties')  

                % Handles for plots
                PDpl_local = zeros(1,3);
                
                subplot(3,1,1)
                PDpl_local(1) = plot(O.calc_time,moments(O.calc_dist,0));                
                ylabel('0^{th} moment [#/g]')
                
                subplot(3,1,2)
                PDpl_local(2) = plot(O.calc_time,moments(O.calc_dist,3));
                ylabel('3^{rd} moment [\mum^3/g]')
                
                subplot(3,1,3)
                PDpl_local(3) = plot(O.calc_time,moments(O.calc_dist,4)./moments(O.calc_dist,3));
                ylabel('Weight average length [\mum]')
                xlabel('Time [s]')
                
                PDpl = [PDpl PDpl_local];
            elseif (~isempty(find(strcmp(plotwhat,'detailed_results'))) || ...
                ~isempty(find(strcmp(plotwhat,'moments'))))
            
                figure(22)
                set(gcf,'numbertitle','off','name','Moments Only')  
                PDpl_local = zeros(1,4);

                subplot(2,2,1)
                PDpl_local(1) = plot(O.calc_time,moments(O.calc_dist,0));
                ylabel('0^{th} moment')
                xlabel('Time')

                subplot(2,2,2)
                PDpl_local(1) = plot(O.calc_time,moments(O.calc_dist,1));
                ylabel('1^{st} moment')
                xlabel('Time')

                subplot(2,2,3)
                PDpl_local(1) = plot(O.calc_time,moments(O.calc_dist,2));
                ylabel('2^{nd} moment')
                xlabel('Time')

                subplot(2,2,4)
                PDpl_local(1) = plot(O.calc_time,moments(O.calc_dist,3));
                ylabel('3^{th} moment')
                xlabel('Time')
                
            end % if
            
            % Process Variables
            if (~isempty(find(strcmp(plotwhat,'results'))) || ...
                ~isempty(find(strcmp(plotwhat,'detailed_results'))) ||...
                ~isempty(find(strcmp(plotwhat,'process'))))
            
                figure(31)
                set(gcf,'numbertitle','off','name','Process Variables (I)')
            
            % Handles for plots
                PDpl_local = zeros(1,1);

                nopvit = 1;
                if ~isempty(O.calc_conc)

                    subplot(2,2,nopvit);                    
                    PDpl_local = plot(O.calc_time,O.calc_conc,'linewidth',1.5);
                    xlim([min(O.calc_time) max(O.calc_time)])
                    xlabel('Time [s]')
                    ylabel('Concentration [g/g]')
                    grid on
                    PDpl = [PDpl PDpl_local];

                    nopvit = nopvit + 1;
                end % if
                
                if ~isempty(O.calc_conc)
                    subplot(2,2,nopvit);                    
                    PDpl_local = plot(O.calc_time,O.calc_conc./O.solubility(O.Tprofile(O.calc_time)),'linewidth',1.5);
                    xlabel('Time [s]')
                    xlim([min(O.calc_time) max(O.calc_time)])
                    ylabel('Supersaturation [-]')
                    grid on
                    PDpl = [PDpl PDpl_local];

                    nopvit = nopvit + 1;
                end % if
                
                if ~isempty(O.calc_conc)
                    subplot(2,2,nopvit);                    
                    PDpl_local = plot(O.calc_time,Tcalc,'linewidth',1.5);
                    xlabel('Time [s]')
                    xlim([min(O.calc_time) max(O.calc_time)])
                    ylabel('Temperature [^\circC]')
                    grid on
                    PDpl = [PDpl PDpl_local];

                    nopvit = nopvit + 1;
                end % if

                if ~isempty(mscalc)

                    subplot(2,2,nopvit);                    
                    PDpl_local = plot(O.calc_time,mscalc,'linewidth',1.5);
                    xlabel('Time')
                    xlim([min(O.calc_time) max(O.calc_time)])
                    ylabel('Total mass Solvent + Antisolvent [g]')
                    grid on
%                     keyboard
                    PDpl = [PDpl PDpl_local(:)'];

                    nopvit = nopvit + 1;
                end

                if ~isempty(O.calc_conc) && ...
                        (~isempty(find(strcmp(plotwhat,'detailed_results'))) || ...
                        ~isempty(find(strcmp(plotwhat,'process'))))
                    PDpl_local = zeros(1,1);

                    figure(32)
                    set(gcf,'numbertitle','off','name','Process Variables (II)')
                    Tvec = linspace(min(Tcalc)-5,max(Tcalc)+5);
                    plot(Tvec,O.solubility(Tvec),'--','linewidth',1.5)
                    legend('Solubility','location','southeast')
                    hold on
                    PDpl_local = plot(Tcalc,O.calc_conc,'linewidth',1.5);
                    xlabel('Temperature')
                    ylabel('Concentration [g/g]')
                    grid on
%                     keyboard
                    PDpl = [PDpl PDpl_local(:)'];
                end
                        
            end % if
            
            if (~isempty(find(strcmp(plotwhat,'detailed_results'))) || ...
                ~isempty(find(strcmp(plotwhat,'integration'))))
                
                if ~isempty(O.calc_conc)
                    figure(41)
                    set(gcf,'numbertitle','off','name',...
                        'Details from Integration')
                    PDpl_local = zeros(1,1);
                    PDpl_local = plot(O.calc_time,massbal(O));
                    xlabel('Time')
                    ylabel('Mass balance [% error]')
                    grid on
                    PDpl = [PDpl PDpl_local];
                end % if
                
                
            end
            
        end % function
        
    end % methods
    
end % classdef