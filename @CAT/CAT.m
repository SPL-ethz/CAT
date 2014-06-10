classdef CAT < hgsetget
    
%% Crystallization Analysis Toolbox (CAT)
%
%% About CAT
%
% The Crystallization Analysis Toolbox is a Matlab toolbox that aims at
% facilitating access to population balance equations.
%
% CAT allows modelling cooling and antisolvent (semi-)batch crystallization
% processes with homogeneous/heterogeneous nucleation and size
% (in)dependent growth, i.e. problems of a general type:
% 
% d(fM)/dt + Md(Gf)/dL = JM*dirac(L)
%
% Where f is a number density function and L a characteristic length of
% the particles; G is a growth rate; M is the total solvent mass. This
% equation is coupled with a mass balance for the solution phase.
%
% CAT is developed at the ETH Zuerich, by Dave Ochsenbein of the Automatic
% Control Laboratory and Martin Iggland of the Separation Processes
% Laboratory.
%
% Website: http://www.ipe.ethz.ch/laboratories/spl/cat
%
%% Usage
%
% CAT uses an object-oriented approach. To use CAT, first define an
% instance of the class by running:
%	>> C = CAT;
%
% CAT objects have several properties which can be set in order to define
% the model you want to solve. To see a list of properties, type
%	>> C
%
% To find out more about any property, use the help function, by typing:
%	>> help CAT.property_name
%
% where property_name is the name of the property you are interested in.
%
% Once all relevant properties have been defined, run the solver:
%	>> C.solve
%
% Analyse results by plotting:
%	>> C.plot
%
% SEE ALSO
% CATTube, Distribution
    
%% CAT Properties

    properties ( Access = protected )
        
        % Time nodes for non-smooth input profiles
        tNodes = [];
        
    end
    
    properties
        
        % The comments before each property are what appear as help, and
        % they also appear in the GUI. Adhere to the following general
        % layout:
        %
        % Property: name
        % Short description of the property
        % Description of how the property is defined. Type, attributes etc,
        % units.
        
        
        % Property: init_dist
        % Initial particle size distribution
        % Defined using Distribution class
        init_dist = Distribution; % Define this default so that C.init_dist.edit works from start!
        
        % Property: init_seed
        % Initial seed mass.
        % Scalar value. The units must be consistent with those used for:
        %  * Initial concentration
        %  * Solubility function
        %  * Crystal density
        init_seed
        
        % Property: init_massmedium
        % Initial mass of the continuous medium: total mass of solvent and
        % antisolvent.
        % Scalar value. The units must be consistent with those used for:
        %  * Initial concentration
        %  * Solubility function
        %  * Antisolvent profile
        %  * Nucleation rate
        init_massmedium
        
        % Property: init_conc
        % Initial concentration.
        % Scalar value, with units of mass solute per total solvent (continuous
        % medium) mass. Use 'S=xx' for a solution with a defined Supersaturation.
        % Units must be consistent with those used for:
        %  * Seed mass
        %  * Crystal density
        %  * Total initial solvent+antisolvent mass
        %  * Solubility function
        %  * Antisolvent profile
        %  * Nucleation rate
        init_conc
        
        
        % Property: Tprofile
        % Temperature profile: temperature as a function of time.
        % Can be: anon. function, matrix or scalar.
        % Units must be consistent with those used for:
        %  * Solubility function
        %  * Nucleation rate
        %  * Growth rate
        %  * Solubility
        Tprofile
        
        % Property: ASprofile
        % Antisolvent profile: mass added as a function of time.
        % Can be: anon. function, or matrix
        % Units must be consistens with those used for:
        %  * Initial concentration
        %  * Solubility function
        %  * Nucleation rate
        %  * Growth rate
        %  * Solvent mass
        ASprofile
        
        
        % Property: solubility
        % Solubility function as a function of temperature (T) and antisolvent mass fraction (xm).
        % Can be an anon. function (Defined as: @(T,xm) or a scalar.
        % Units must be mass solute per total solvent mass and consistent with:
        %  * Seed mass
        %  * Crystal density
        %  * Total initial solvent+antisolvent mass
        %  * Initial concentration
        %  * Antisolvent profile
        %  * Nucleation rate
        solubility
        
        % Property: rhoc
        % Crystal density.
        % Scalar value. The units must be consistent with those used for:
        %  * Seed mass
        %  * Initial concentration
        %  * Solubility function
        rhoc
        
        % Property: kv
        % Shape factor
        % Dimensionless scalar value, between 0 and 1.
        kv
        
        % Property: growthrate
        % Growth rate as a function of supersaturation (S), temperature (T) and size (y).
        % Defined as an anon. function: @(S,T,y)
        % Units must be consistent with those used for:
        %  * Initial distribution
        %  * Temperature profile
        %  * Antisolvent profile
        %  * Nucleation rate
        %  * Solution time
        %
        % The function should return a vector the same size as y
        growthrate
        
        % Property: nucleationrate
        % Nucleation rate function
        % Nucleation rate as a function of supersaturation (S), temperature (T) or (moments of) the distribution F.
        % Defined as an anon. function: @(S,T,F)
        % Units must be consistent with those used for:
        %  * Initial distribution
        %  * Temperature profile
        %  * Antisolvent profile
        %  * Growth rate
        %  * Solution time
        %
        % The function should return a scalar value
        nucleationrate
        
        
        % Property: sol_time
        % Solution time
        % Vector or scalar.
        % The units must be consistent with those used for:
        %  * Growth rate
        %  * Nucleation rate
        %  * Antisolvent profile
        %  * Temperature profile
        sol_time
        
        % Property: sol_method
        % Solution method.
        % Defines which numerical method to use.
        % Use the solutionMethods function to see a list of available
        % solvers.
        sol_method
        
        % Property: sol_options
        % Solver options
        sol_options = [];
        
        %
        % Results properties
        %
        
        % Property: calc_time
        % Vector of actual times returned by solver
        calc_time
        
        % Property: calc_dist
        % Distributions for each time step
        calc_dist
        
        % Property: calc_conc
        % Vector of solution concentrations for times given in calc_time
        calc_conc
  
    end % properties
    
    %% CAT Methods
    
    methods

        %% - Set and Get methods

        %% -- set.rhoc
        
        function set.rhoc(O,value)
            
            % SET.rhoc
            %
            % Check the crystal density. It must be a scalar
            
            if isempty(value) || O.diagnose('rhoc',value)
            
                O.rhoc = value;
                
            end % if else
            
            O.rhoc_onset;
            
        end % function
        
        %% -- set.kv
        
        function set.kv(O,value)
            
            % SET.kv
            %
            % Check the shape factor. It must be a scalar
            
            if isempty(value) || O.diagnose('kv',value) 
            
                O.kv = value;

            end % if else
            
            O.kv_onset;
            
        end % function
        
        %% -- set.init_dist
        
        function set.init_dist(O,value)
            
            % SET.INIT_DIST
            %
            % Check the initial distribution, it must be a distribution
            % class object
            
            if isempty(value) || O.diagnose('init_dist',value) 
            
                O.init_dist = value;

                
            end % if else
            
            
            % Extra function - overwritable in subclasses
            O.init_dist_onset;
            
        end % function
        
        %% -- set.init_seed
        
        function set.init_seed(O,value)
            
            % SET.INIT_SEED
            %
            % Set mass of seeds
            
            if isempty(value) || O.diagnose('init_seed',value)
                
                O.init_seed = value;
                
            end
            
            % Extra function - overwritable in subclasses
                O.init_seed_onset;
            
        end % function
        
        %% -- set/get.init_conc
        
        function set.init_conc(O,value)
            
            % SET.INIT_CONC
            %
            % Check the initial concentration, it must be a positive,
            % finite scalar (can be zero)
            if isempty(value) || O.diagnose('init_conc',value)
                O.init_conc = value;
                
            end
            
            % Extra function - overwritable in subclasses
            O.init_conc_onset;
            
        end % function
        
        function [cinit] = get.init_conc(O)
            
            % GET.INIT_CONC
            %
            % Getter method for init conc

            if ischar(O.init_conc) && ~isempty(strfind(O.init_conc,'S=')) && ~isempty(O.solubility)
                
                S0 = str2double(strrep(O.init_conc,'S=',''));
                if isnan(S0)
                    S0 = eval(strrep(O.init_conc,'S=','')); % maybe user has written something 'S=2/3'
                end
                if ~isempty(O.solubility) && ~isempty(O.sol_time(1)) && ~isempty(O.Tprofile) && ~isempty(O.ASprofile) && ~isempty(O.init_massmedium) && isa(O.solubility,'function_handle')
                    O.init_conc = evalanonfunc(O.solubility,O.Tprofile(O.sol_time(1)),O.ASprofile(O.sol_time(1))/O.init_massmedium)*S0;
                end
            end
            cinit = O.init_conc;
            
        end % function
        
        %% -- set.init_massmedium
        
        function set.init_massmedium(O,value)
            
            if isempty(value) || O.diagnose('init_massmedium',value)
                
                O.init_massmedium = value;

            end
            
            % Extra function - overwritable in subclasses
                O.init_massmedium_onset;
            
        end % function
        
        %% -- set.solubility
        
        function set.solubility(O,value)
            
            % SET.solubility
            %
            %  Setter method for solubility Must be a function handle with
            %  1 or 2 inputs
            
%             % check if the value can be easily calculated to a number (a
%             % number may come in form of a string from the GUI)
%             try %#ok<TRYNC>
%                 value = str2num(value); %#ok<*ST2NM>
%             end
            
            % Check for number - convert to constant function
            if isempty(value) || O.diagnose('solubility',value)
                if isnumeric(value) && length(value) == 1
                    O.solubility = str2func(['@(T,xm)' num2str(value) '*ones(size(T))']);
                elseif ischar(value)
                    % Check for string 
                    if isempty(strfind(value,'@'))
                        O.solubility = str2func(['@(T,xm)' value '*ones(size(T))']);
                    else
                        O.solubility = str2func([value '*ones(size(T))']);
                    end
                elseif isa(value,'function_handle') || isempty(value)
% 
%                     if isa(value,'function_handle') && nargin(value) < 2
%                         % This is too few, the function needs to accept two
%                         % inputs (even if the second isn't used)
%                         O.solubility = str2func(['@(T,xm)' anonfunc2str(value)]);
%                     else
%                         % If the function is defined with more inputs, there is
%                         % no problem as long as the other values are not
%                         % needed. This will give an error later
                        O.solubility = value;
%                     end % if elseif
%                 
                end % if else
% 
            end
            
            % Extra function - overwritable in subclasses
                O.solubility_onset;
            
        end % function
        
        %% -- set.sol_time
        
        function set.sol_time(O,value)
            
            % SET.SOL_TIME
            %
            % Check the solution time vector. This should be a vector of
            % monotonically increasing values
            
            if isempty(value) || O.diagnose('sol_time',value)
            
                if isempty(value) || length(value) > 1 
                    O.sol_time = value(:)';
                elseif isscalar(value)
                    O.sol_time = [0 value];
                end % if else
                
            end
            
            
            % Extra function - overwritable in subclasses
                O.sol_time_onset;
        end % function
        
        %% -- set.sol_method
        
        function set.sol_method(O,value)
            
            % SET.SOL_METHOD
            %
            % Check the defined solution method. Should be a string
            %
            % Probable more checks should be carried out at this point in
            % the future
            
            if O.diagnose('sol_method',value)
                if ischar(value)

                    % Remove spaces from name first
                    value = strrep(value,' ','');

                    % don't be case sensitive and allow alternative forms
                    switch lower(value)
                        case {'cd','centraldifference'}
                            O.sol_method = 'centraldifference';
                        case {'mp','movingpivot'}
                            O.sol_method = 'movingpivot';
                        case {'hr','hires','highresolution'}
                            O.sol_method = 'hires';
                        otherwise
                            error('CAT:SetSol_Method:unknown',...
                                'Unknown solution method');
                    end % switch

                elseif isempty(value)
                    warning('CAT:SetSol_Method:isempty',...
                        'The property sol_method was set to the default value (centraldifference)');
                    O.sol_method = 'centraldifference';

                end % if else
            
            
            end
            
            % Extra function - overwritable in subclasses
                O.sol_method_onset;
            
        end % function
        
        %% -- set.sol_options
        
        function set.sol_options(O,value)
            
            % SET.SOL_OPTIONS
            %
            % Check the defined OPTIONS. Should be a cell array
            %
            % Probable more checks should be carried out at this point in
            % the future
            
            if O.diagnose('sol_options',value)
                O.sol_options = value;

             
            end
            
            % Extra function - overwritable in subclasses
                O.sol_options_onset;
            
        end % function
        
        
        %% -- set.Tprofile
        
        function set.Tprofile(O,value)
            
            % SET.Tprofile
            %
            % Check the Temperature profile. It must be a matrix with
            % positive, finite elements. The first row indicates the times
            % of the nodes whereas the second row indicates Temp's

            % check if the value can be easily calculated to a number (a
            % number may come in form of a string from the GUI)
            try %#ok<TRYNC>
                value = str2num(value);
            end
            
            if isempty(value) || O.diagnose('Tprofile',value)
                if ~isempty(value) && ismatrix(value) && length(value(:,1))==2 && all(isfinite(value(:)))

                    if isempty(O.sol_time)
                        O.sol_time = value(1,end);
                    end
                    
                    if value(1,end)<O.sol_time(end)
                        value = [value [0;0]];
                        value(1,end) = O.sol_time(end);
                        value(2,end) = value(2,end-1);
                    end
    
                    O.Tprofile = str2func(strcat('@(t) piecewiseLinear(',data2str(value(1,:)),',',data2str(value(2,:)),',t)')); %
                    O.tNodes = unique([O.tNodes value(1,:)]);
                
                elseif isempty(value)
                    O.Tprofile = value;

                elseif isa(value,'function_handle') && nargin(value)==1

                    if length(value([0 10]))==1
                        O.Tprofile = @(t) value(t)*ones(size(t)); % output should be a vector of same size as input
                    else
                        O.Tprofile = value;
                    end

                elseif isnumeric(value) && isscalar(value)
                    O.Tprofile =  str2func(strcat('@(t)', data2str(value),'*ones(size(t))'));

                

                end % if else

            
            end
            
            % Extra function - overwritable in subclasses
                O.Tprofile_onset;
            
        end % function
        
        %% -- set.ASprofile
        
        function set.ASprofile(O,value)
            
            % SET.ASprofile
            %
            % Check the time profile (for added AS profiles). It must be
            % a strictly increasing(!), positive vector
            
            % check if the value can be easily calculated to a number (a
            % number may come in form of a string from the GUI)
            try %#ok<TRYNC>
                value = str2num(value);
            end
            
            if isempty(value) || O.diagnose('ASprofile',value)
                if ~isempty(value) &&  ismatrix(value) && length(value(:,1))==2 && all(isfinite(value(:))) && all(diff(value(2,:))>=0)

                    if value(1,end)<O.sol_time(end)
                        value = [value [0;0]];
                        value(1,end) = O.sol_time(end);
                        value(2,end) = value(2,end-1);
                    end
                    O.ASprofile = @(t) piecewiseLinear(value(1,:),value(2,:),t); %
                    O.tNodes = unique([O.tNodes value(1,:)]);
                
                elseif isempty(value)
                    O.ASprofile = value;

                elseif isa(value,'function_handle') && nargin(value)==1
                    O.ASprofile = value;
                    
                elseif isnumeric(value) && isscalar(value)
                    O.ASprofile =  str2func(strcat('@(t)', data2str(value),'*ones(size(t))'));

                elseif isempty(value)
                    O.ASprofile = value;
                else
                    warning('CAT:SetTprofile:WrongType',...
                        'The ASprofile property must be a positive, finite matrix (may be zero) or a function handle with one input');

                end % if else

                
            end
            
            % Extra function - overwritable in subclasses
                O.ASprofile_onset;
            
        end % function
        
        %% -- set.tNodes
        
        function set.tNodes(O,value)
            
            % SET.tNodes
            %
            % Set time nodes (make sure integrator covers them correctly)
            if isvector(value) && all(value>=0) || isempty(value)

                O.tNodes = value;

            else
                warning('CAT:SettNodes:WrongType',...
                    'The tNodes property must be a non-negative vector');
                
            end % if else
            
        end % function
        
        %% -- set.growthrate
        
        function set.growthrate(O,value)
            
            % SET.GROWTHRATE
            %
            % Check the growth rate: should be a function handle, accept 3
            % arguments: S (scalar), T (temperature), and y (vector). The output should be
            % the same size as y
            
            % If the growth rate is not given as a function, make best
            % choice to convert it into 
            
            % check if the value can be easily calculated to a number (a
            % number may come in form of a string from the GUI)
            try %#ok<TRYNC>
                value = str2num(value);
            end
            
            if isempty(value) || O.diagnose('growthrate',value)
                if isnumeric(value) && length(value) == 1
                    O.growthrate = str2func(['@(S,T,y)' num2str(value) '*ones(size(y))']);
                elseif isempty(value)
                    O.growthrate = value;
                elseif ischar(value)
                    if isempty(strfind(value,'@'))
                        O.growthrate = str2func(['@(S,T,y)' value]);
                    else
                        O.growthrate = str2func(value);
                    end
                elseif isa(value,'function_handle')

                    % Check the number of inputs
                    if nargin(value) == 3

                        % Check the output using 2 example values
                        out = value(1.1,1,linspace(0.1,1,10));

                        % Check size
                        if any( size(out) ~= [1 10] )
                            % Size of output wrong
                            warning('Distribution:setgrowthrate:Wrongsize',...
                                'The growth rate function returns a vector which is not the same size as the input vector');
                        end
                        % Set the growthrate anyway
                        O.growthrate = value;

                    elseif nargin(value) == 2 && length(value(1.1,1))==1
                        % assume size independent function
                        
                        y = strsplit(data2str(value),')');
                        f = [y{1},',y)',implode(y(2:end),')'),'*ones(size(y))'];
                        O.growthrate = str2func(f);

                    elseif nargin(value) == 2 && length(value(1.1,[1 2]))==2
                        % assume temperature independent function
                        O.growthrate = @(S,~,y) value(S,y);

                    elseif nargin(value) == 1
                        % assume growth rate function only depending on S
                        y = strsplit(data2str(value),')');
                        f = [y{1},',~,y)',implode(y(2:end),')'),'*ones(size(y))'];
                        O.growthrate = str2func(f);
                    else
                        warning('Distribution:setgrowthrate:Wrongnargin',...
                            'The growth rate function must have 3 input arguments (supersaturation, temperature, sizes) or 2 input arguments (S,T) or (S,y)');
                    end % if
                
                end %if
                
            end
            
            % Extra function - overwritable in subclasses
                O.growthrate_onset;
            
        end % function
        

        %% -- set.nucleationrate
        
        function set.nucleationrate(O,value)
            
            % SET.nucleationrate
            %
            % Check the nucleationrate rate: should be a function handle,
            % accept max. 3 arguments: S (supersaturation), T (temperature), F (distribution). The output should be
            % a scalar
            if O.diagnose('nucleationrate',value)
                if isnumeric(value) && length(value) == 1
                    O.nucleationrate = str2func(['@(S,T,F)' num2str(value) '*ones(size(S))']);
                elseif ischar(value)
                    if isempty(strfind(value,'@')) % check whether string is already in complete an. function form
                        O.nucleationrate = str2func(['@(S,T,F)' value]);
                    else
                        O.nucleationrate = str2func(value);
                    end

                elseif isa(value,'function_handle')
                    if nargin(value) == 1
                        O.nucleationrate = str2func(['@(S,~,~)' anonfunc2str(value)]);
                    elseif nargin(value) == 2
                        O.nucleationrate = str2func(['@(S,T,~)' anonfunc2str(value)]);
                    else
                        O.nucleationrate = value;
                    end

                elseif isempty(value)
                    O.nucleationrate = value;

                end %if

            
            end
            
            % Extra function - overwritable in subclasses
                O.nucleationrate_onset;
            
        end % function

        
        %% - Additional CAT methods
        
        %% -- massmedium: Mass solvent + antisolvent at t
        function mscalc = massmedium(O,t)
            if ~exist('t','var')
                t = O.calc_time(:);
            end
            mscalc = O.init_massmedium+O.ASprofile(t)-O.ASprofile(0); % total amount of medium  
        end % function
        
        %% -- massbal
        
        function PDma = massbal(O)
            mass_solute = O.calc_conc(:).*massmedium(O);   % total mass of solute
            m3 = moments(O.calc_dist,3);                            % third moment over time
            mass_crystals = O.rhoc*O.kv*massmedium(O).*m3(:);    % total mass of crystals
            
            PDma = 100*((mass_solute + mass_crystals)/(mass_solute(1)+mass_crystals(1))-1);     % mass balance error   
        end % function
        
    end % methods
    
    %% Hidden methods
    
    %% - All onset methods
    % Do nothing - no function in this class, merely something to
    % be overwritten by subclass
    methods (Hidden)
        
        function init_dist_onset(O) %#ok<*MANU>
            
        end % function
        
        
        function init_conc_onset(O)
            
        end % function
        
        
        function init_seed_onset(O)
            
        end % function
        
        
        function init_massmedium_onset(O)
            
        end % function
        
        
        function rhoc_onset(O)
            
        end % function
        
        
        function kv_onset(O)
            
        end % function
        
        
        function solubility_onset(O)
            
        end % function
        
        
        function Tprofile_onset(O)
            
        end % function
        
        
        function ASprofile_onset(O)

        end % function
        
        
        function growthrate_onset(O)
            
        end % function
        
        
        function nucleationrate_onset(O)

        end % function
        
        
        function sol_time_onset(O)
            
        end % function
        
        
        function sol_method_onset(O)
            
        end % function
        
        
        function sol_options_onset(O)
            
        end % function
    end
    
    %% Static Methods
    methods (Static)
        % Compare two CAT objects
        [results,same,different] = compare(CAT1,CAT2,varargin); 
            
        % Fill out the template form
        fillOutForm(setupCat);
        
    end
    
end % classdef