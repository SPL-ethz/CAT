classdef CAT < hgsetget
    %% Class CAT
    % The Crystallization Analysis Toolbox (CAT) class defines numerical
    % settings, initial and operating conditions to set up problems
    % encountered in crystallization science. It furthermore allows to
    % solve this kind of problems using a variety of different solvers.
    % Construct a CAT object and set properties as you see fit. The
    % currently recommended pathway is to fill out protocat and save it as
    % new file.
    
    % Dave Ochsenbein, 15.10.2013
    
    properties ( Access = protected )
        
        % nodes for non-smooth input profiles
        tNodes = [];
        
    end
    
    properties
        
        % Initial distribution (Distribution object, defines size and
        % distribution
        init_dist
        
        % Initial concentration
        init_conc % saturated
        
        % Solubility
        % @(T,xm)
        solubility
        
        % Temperature profile
        Tprofile
        
        % (Anti)solvent added mass profile
        ASprofile

        % Growth rate function This function should be called as
        % growthrate(S,T,y) where S is the current supersaturation, T is the
        % temperature and y is the size. It should return a vector the same
        % size as y
        growthrate
        
        % Nucleation rate function
        % This function should be called as nucleationrate(c,T,m) where S is the current
        % supersaturation, T is the temperature. Optionally, the user can
        % specificy that the nucleation rate depends on a moment m of the
        % passed distribution F
        nucleationrate
        
        % Seed mass
        init_seed
        
        % Initial mass of solvent + antisolvent
        init_massmedium
        
        % Solution time vector
        sol_time
        
        % Crystal density
        rhoc
        
        % Shape factor
        kv
        
        % Method to use - default to central difference
        sol_method
        
        % Solver Options
        sol_options = [];
        
        %% Results
        
        % Vector of actual times returned by solver
        calc_time
        
        % Distributions for each time step
        calc_dist
        
        % Concentrations over time [g solute / g total solvents]
        calc_conc
  
    end % properties
    
    methods
        
        %% Method CAT (constructor)
        
        function O = CAT(varargin)
            
            % CAT
            %
            % Constructor function for CAT class. Called as:
            %	PD = CAT(init_dist,init_conc,sol_time,...
            %      sol_method,growthrate)
            %
            % Leave any of the variables empty to not assign them.
            %
            % The variables which are given are assigned directly to the
            % corresponding property. All of them can also be assigned
            % later.
            
            % Set values, if they are given
            
            if  nargin == 0 || (nargin>0 && ~isempty(find(strcmp(varargin,'empty'), 1)))
                
            elseif  (nargin>0 && ~isempty(find(strcmp(varargin,'default'), 1)))
                
                O.setDefaults;
                
            end
            
        end % function        
        
%% SETTER AND GETTER METHODS

        %% Method set.rhoc
        
        function set.rhoc(O,value)
            
            % SET.rhoc
            %
            % Check the crystal density. It must be a scalar
            
            if isempty(value) || O.diagnose('rhoc',value)
            
                O.rhoc = value;
                
            end % if else
            
            O.rhoc_onset;
            
        end % function
        
        %% Method set kv
        
        function set.kv(O,value)
            
            % SET.kv
            %
            % Check the shape factor. It must be a scalar
            
            if isempty(value) || O.diagnose('kv',value) 
            
                O.kv = value;

            end % if else
            
            O.kv_onset;
            
        end % function
        
        %% Method set init_dist
        
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
        
        %% Method set.init_seed
        
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
        
        %% Method set/get init_conc
        
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
            
            
            if strcmpi('sat',O.init_conc) && ~isempty(O.solubility)
                O.init_conc = O.solubility(O.Tprofile(O.sol_time(1)),O.ASprofile(O.sol_time(1))/O.init_massmedium);
            end
            cinit = O.init_conc;
            
        end % function
        
        %% Method set.init_massmedium
        
        function set.init_massmedium(O,value)
            
            if isempty(value) || O.diagnose('init_massmedium',value)
                
                O.init_massmedium = value;

            end
            
            % Extra function - overwritable in subclasses
                O.init_massmedium_onset;
            
        end % function
        
        %% Method set.solubility
        
        function set.solubility(O,value)
            
            % SET.solubility
            %
            %  Setter method for solubility Must be a function handle with
            %  1 or 2 inputs
            
            % check if the value can be easily calculated to a number (a
            % number may come in form of a string from the GUI)
            try %#ok<TRYNC>
                value = str2num(value); %#ok<*ST2NM>
            end
            
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

                    if isa(value,'function_handle') && nargin(value) < 2
                        % This is too few, the function needs to accept two
                        % inputs (even if the second isn't used)
                        O.solubility = str2func(['@(T,xm)' anonfunc2str(value)]);
                    else
                        % If the function is defined with more inputs, there is
                        % no problem as long as the other values are not
                        % needed. This will give an error later
                        O.solubility = value;
                    end % if elseif
                
                end % if else

            end
            
            % Extra function - overwritable in subclasses
                O.solubility_onset;
            
        end % function
        
        %% Method set.sol_time
        
        function set.sol_time(O,value)
            
            % SET.SOL_TIME
            %
            % Check the solution time vector. This should be a vector of
            % monotonically increasing values
            
            if isempty(value) || O.diagnose('sol_time',value)
            
                if isempty(value) || length(value) > 1 
                    O.sol_time = value;
                elseif isscalar(value)
                    O.sol_time = [0 value];
                end % if else
                
            end
            
            
            % Extra function - overwritable in subclasses
                O.sol_time_onset;
        end % function
        
        %% Method set.sol_method
        
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
        
        %% Method set.sol_options
        
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
        
        
        %% Method set.Tprofile
        
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
        
        %% Method set.ASprofile
        
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
    %                 O.ASprofile = @(t) interp1(value(1,:),value(2,:),t); %
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
        
        %% Method set.tNodes
        
        function set.tNodes(O,value)
            
            % SET.tNodes
            %
            % Set time nodes (make sure integrator covers them correctly)
            if isvector(value) && all(value>=0)

                O.tNodes = value;

            else
                warning('CAT:SettNodes:WrongType',...
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
                        O.growthrate = @(S,T,~) value(S,T);

                    elseif nargin(value) == 2 && length(value(1.1,[1 2]))==2
                        % assume temperature independent function
                        O.growthrate = @(S,~,y) value(S,y);

                    else
                        warning('Distribution:setgrowthrate:Wrongnargin',...
                            'The growth rate function must have 3 input arguments (supersaturation, temperature, sizes) or 2 input arguments (S,T) or (S,y)');
                    end % if
                
                end %if
                
            end
            
            % Extra function - overwritable in subclasses
                O.growthrate_onset;
            
        end % function
        

        %% Method set.nucleationrate
        
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

        
%% Additional CAT methods
        %% Mass solvent + antisolvent at t
        function mscalc = massmedium(O,t)
            if ~exist('t','var')
                t = O.calc_time(:);
            end
            mscalc = O.init_massmedium+O.ASprofile(t)-O.ASprofile(0); % total amount of medium  
        end % function
        
        %% Method massbal
        
        function PDma = massbal(O)
            mass_solute = O.calc_conc(:).*massmedium(O);   % total mass of solute
            m3 = moments(O.calc_dist,3);                            % third moment over time
            mass_crystals = O.rhoc*O.kv*massmedium(O).*m3(:);    % total mass of crystals
            
            PDma = 100*((mass_solute + mass_crystals)/(mass_solute(1)+mass_crystals(1))-1);     % mass balance error   
        end % function
        
        %% Method save
        % This function saves the CAT instance in a mat file
        function save(O,CATname)
            
            if nargin<2 || isempty(CATname{1})
                namestr = 'kitten';
            else
                namestr = CATname{1};
            end
            
            namestr = strcat(namestr,'_',datestr(now,'YYYYmmdd_HHMMSS'),'.mat');
 
            superCAT.kitty = O; %#ok<STRNU>
            save(namestr,'-struct','superCAT')
            disp(strcat({'Saved file to '},namestr))
            clear superCAT n
        end % function
        
            
        %% Method clone
        % This function clones the CAT instance
        function [copyCAT] = clone(O,Original)
            
            
            if ~isa(Original,'CAT') && ~isa(Original,'CATTube')
               
                warning('CAT:clone:notaCAT',...
                    'The object you try to clone must be of class CAT');
                
            else
                
                if nargout == 0
                    F = O;
                else
                    if isa(Original,'CATTube')
                        F = CATTube('hush','uncloned');
                    elseif isa(Original,'CAT')
                        F = CAT('hush');
                    end
                end
                    
                
                fieldnames = properties(O);
                
                for i = 1:length(fieldnames)
                   
                    F.(fieldnames{i}) = Original.(fieldnames{i});
                    
                end
                    
                if nargout == 1
                    copyCAT = F;
                else
                    copyCAT = [];
                end
                
            end
            
            
        end % function
        
        
        
    end % methods
    
    %% All onset methods
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