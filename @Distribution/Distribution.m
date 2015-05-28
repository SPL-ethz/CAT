classdef Distribution < Easyset & matlab.mixin.Copyable
    %% Class Distribution
    %
    %% About Distribution class
    %
    % The Distribution class defines particle size distributions as a
    % function of a characteristic variable (normally a length). It handles
    % both analytic and numerical distributions.
    %
    % The primary use of this class is within the Crystallisation Analysis
    % Toolbox (CAT).
    %
    %% Usage
    %
    % Note that the distribution is assumed to be a density function and
    % therefore boundaries of the bins are necessary to correctly calculate
    % the moments. If boundaries are not assigned, it will be assumed that
    % boundary(1) = 0 and that the remaining pivot sizes (y) lie in the
    % arithmetic mean of the remaining boundary points.
    %
    % The size vector needs to be defined. The distribution can be defined
    % as either a function handle or as a vector. The class checks the
    % inputs, and returns a vector in any case.
    %
    % The class can be called as follows:
    %
    %   dist = Distribution(y,@(x)normalpdf(x,1,0.1));
    %   or
    %   dist = Distribution(y,{'normal',1,0.1});
    %
    % to create a distribution based on the normal distribution, or the
    % properties y and F can be set separately:
    %
    %   dist = Distribution;
    %   dist.y = linspace(0,2,50);
    %   dist.F = @(x) normpdf(x,1,0.1);
    %
    % This class can be used to define multiple distributions, having their
    % own size vectors, in the following way:
    %
    %   dist(1) = Distribution(linspace(0,2),@(x)normpdf(x,1,0.1));
    %   dist(2) = Distribution(linspace(0,3),@(x)normpdf(x,2,0.5));
    %
    % Moments of Distribution objects can be calculated using the moments
    % method. See the corresponding entry.
    %
    % SEE ALSO
    % CAT, CATTube, Distribution.moments, Distribution.plot,
    % Distribution.dist2str, Distribution.data2str
    
    %% Properties
    
    % The comments before each property are what appear as help, and
    % they also appear in the GUI. Adhere to the following general
    % layout:
    %
    % Property: name
    % Short description of the property
    % Description of how the property is defined. Type, attributes etc,
    % units.
    
    properties
        
        % Property: y
        % Vector of coordinates at which distribution is defined (pivots)
        % Real, positive vector
        y = linspace(0,1,11);
        
        % Property: boundaries
        % Vector of coordinates at which distribution is defined (boundaries)
        % Real, positive vector
        boundaries
        
    end % properties
    
    properties (Dependent)
        
        % Property: F
        % Values of the particle size distribution
        % Either an anonymous function or a vector
        F = [];
        
    end % properties
    
    properties (Access=private)
        
        % Property: pF
        % Private version of F - this is where the data is stored.
        pF = [];
        
        % Property: mu
        % Mean of the distribution (if it is defined in such a
        % way).
        mu = [];
        
        % Property: sigma
        % Mean of the distribution (if it is defined in such a
        % way).
        sigma = [];
        
    end % properties
    
    %% Methods
    
    methods
        
        %% - Distribution (constructor)
        
        function O = Distribution(y,F,boundaries)
            
            % Distribution constructor
            %
            % If given, define y, F and boundaries based on inputs.
            % Otherwise, do nothing
            %
            % SEE ALSO
            % Distribution
            
            if nargin > 0 && ~isempty(y)
                O.y = y;
            end % if
            
            if nargin > 1 && ~isempty(F)
                O.F = F;
            end % if
            
            if nargin > 2 && ~isempty(boundaries)
                O.boundaries = boundaries;
            elseif nargin == 2 || (nargin==3 && isempty(boundaries))
                O.boundaries = [0 (O.y(2:end)+O.y(1:end-1))/2 3/2*O.y(end)-O.y(end-1)/2];
            end % if
            
        end % function
        
        %% Method set.y,
        
        % Call the checkPropertyValue() method
        % Along with the set() method, this function
        % captures the two ways of setting a property:
        % D.y = ...
        % set(D,'y',...)
        
        function set.y(O,value)
            
            % Distribution.set.Y
            %
            % Check input for y: a vector, all positive, increasing. Always
            % set y as a row vector
            %
            % SEE ALSO
            % Distribution
            
            % Define properties here
            O.classes.y = 'numeric';
            O.attributes.y = {'vector','real','finite','nondecreasing'};
            
            % Redirect to checkPropertyValue function to do checking
            O.checkPropertyValue('y',value);
            
            % If checks didn't fail, method has not exited, so do extra
            % checks
            

            O.y = value(:)';
                

            
        end % function
        
        %% Method set.boundaries
        
        % Call the checkPropertyValue() method
        % Along with the set() method, this function
        % captures the two ways of setting a property:
        % D.y = ...
        % set(D,'y',...)
        
        function set.boundaries(O,value)
            
            % Distribution.set.boundaries
            %
            % Check input for boundaries: a vector, all positive (and 0),
            % duplicates allowed. Always set boundaries as a row vector
            %
            % SEE ALSO
            % Distribution
            
            % Define properties here
            O.classes.boundaries = 'numeric';
            O.attributes.boundaries = {'vector','real','finite','nondecreasing'}; 
            
            % Redirect to checkPropertyValue function to do checking
            O.checkPropertyValue('boundaries',value);
            
            % If checks didn't fail, method has not exited, so do extra
            % checks
            
            if length(value) > 1 
                O.boundaries = value(:)';
                
                if length(O.y) ~= length(value)-1
                    O.y = (O.boundaries(1:end-1)+O.boundaries(2:end))/2;
%                     warning('Distribution:setboundaries:yAndBoundariesConsistency',...
%                         'Numel of property boundaries inconsistent with property y. Reset y as arithmetic means of boundaries.');
                end

            else
                warning('Distribution:setboundaries:yIsScalar',...
                    'Property boundaries must be a vector with at least two elements');

            end % if
            
        end % function
        
        %% - set.F
        
        function set.F(O,value)
            
            % Distribution.set.F
            %
            % Check input: vector, function handle, empty
            %
            % SEE ALSO
            % Distribution
            
            % Define properties here
            O.classes.F = {'function_handle','numeric','cell'};
            O.attributes.F = {'vector','real'};
            
            % Redirect to checkPropertyValue function to do checking
            O.checkPropertyValue('F',value);
            
            if isa(value,'function_handle')
                % Value is ok, set
                O.pF = value;
            elseif ~iscell(value) && isvector(value) || isempty(value)
                O.pF = value(:)'; %make it a row vector
            elseif iscell(value) && length(value)==3
                if strcmpi(value{1},'normal')
                    O.pF = str2func(['@(x) 1./(',data2str(value{3}),'*sqrt(2*pi))*exp(-((x-',data2str(value{2}),').^2/(2*',data2str(value{3}),'^2)))']);
                    O.mu = value{2};
                    O.sigma = value{3};
                elseif strcmpi(value{1},'lognormal')
                    O.pF = str2func(['@(x) 1./(x*',data2str(value{3}),').*exp(-((log(x)-',data2str(value{2}),').^2/(2*',data2str(value{3}),'^2)))']);
                    O.mu = value{2};
                    O.sigma = value{3};
                end
            else
                % Value is not OK - display warning
                warning('Distribution:setF0:WrongType',...
                    'F has to be a function handle, a vector or a cell of length 3');
            end % if
            
        end % function
        
        %% - get.F
        
        function Fout = get.F(O)
            
            % Distribution.get.F
            % Return F, based on what is in y.
            % If F is a function handle: calculate vector using positions
            % at y. If F is a vector, check so that its size is the same
            % as that of y
            %
            % SEE ALSO
            % Distbribution
            
            if isa(O.pF,'function_handle')
                % Function handle - calculate values
                
                Fout = O.pF(O.y);
                
            elseif isvector(O.pF)
                
                if any( size(O.pF) ~= size(O.y) )
                    warning('Distribution:getF0:WrongSize',...
                        'F is not the same size as y');
                    Fout = O.pF;
                else
                    Fout = O.pF;
                end % if else
                
            else
                Fout = [];
            end % if else
            
            % Always return a row vector
            Fout = Fout(:)';
            
        end % function
        
        %% - getFunction
        
        function fnc = getFunction(O)
            
            % Distribution.getFunction
            %
            % Return the actual function definition of the Distribution -
            % not just the values. If no function is defined, return empty
            %
            % SEE ALSO
            % Distribution
            
            if isa(O.pF,'function_handle')
                fnc = O.pF;
            else
                fnc = [];
            end % if
            
        end % function
        
        %% - moments(F,j)
        
        function Fmo = moments(O,j,icalc)
            
            % Distribution.moments
            %
            % Calculate moments from distribution/distributions over time.
            %
            % Use moments(F,j) to calculate the j'th moment of F over its
            % length.
            %
            % moments returns the indicated moment as a row vector
            %
            % To calculate moments at specific indices points, use e.g.:
            %   F.moments(3,[1 3 5])
            %
            % SEE ALSO
            % Distribution
            
            
            if nargin <= 2 && ~exist('icalc','var')
                icalc = 1:length(O);
            end
            
            % If possible, calculate moment. Otherwise do nothing
            if nargin > 1 && ~isempty(j)
                
                % Preallocate Fmo
                Fmo = zeros(size(icalc));
                
                for i = icalc
                    
                    if ~isempty(O(i).F)
                        
                        if isempty(O(1).boundaries)
                            Dy = diff([0 O(i).y]);
                        else
                            Dy = diff(O(i).boundaries);
                        end
                        
                        Fmo(icalc==i) = sum(O(i).F(:) .* Dy(:).* O(i).y(:).^j);
                        
                    end % if
                    
                end % for
                
            else
                warning('Distribution:moments:nomoment',...
                    'No type of moment was indicated');
                Fmo = [];
            end % end
            
            % Always return a row vector
            Fmo = Fmo(:)';
            
        end % function
        
        %% - disp
        
        function disp(O)
            
            % Distribution.disp
            %
            % Display the distribution in a string representation
            %
            % SEE ALSO
            % Distribution
            
            fprintf([dist2str(O) 10]);
            
        end % function
        
        %% - data2str
        
        function outstr = data2str(O)
            
            % Distribution.data2str
            %
            % Returns the distribution as a string (useful for a comparison of
            % distributions, which in general can be vectors or function
            % handles).
            %
            % This string can be used to recreate the distribution
            %
            % SEE ALSO
            % Distribution, Distribution.dist2str
            
            if isnumeric(O.pF)
                Fstr = mat2str(O.pF);
            elseif isa(O.pF,'function_handle')
                Fstr = func2str(O.pF);
                if ~isempty(O.sigma)
                    Fstr = strrep(Fstr,'value{2}',num2str(O.mu));
                    Fstr = strrep(Fstr,'value{3}',num2str(O.sigma));
                end
            end
            outstr = strcat('Distribution(',data2str(O.y),',',Fstr,',',data2str(O.boundaries),')');
            
        end % function data2str
        
        %% - dist2str
        
        function outstr = dist2str(O)
            
            % Distribution.dist2str
            %
            % Returns a string representation of the distribution object,
            % with shorthand notation describing main properties
            %
            % This string can not be used to recreate the distribution -
            % use data2str for this
            %
            % SEE ALSO
            % Distribution, Distribution.data2str
            
            % Returns a string representation of the distribution
            if isa(O.pF,'function_handle')
                type = 'Fnc';
                outstr = sprintf('%s; d_10 = %.2g, m_3 = %.2g',...
                    type,O.moments(1)/O.moments(0),O.moments(3) );
            elseif isvector(O.pF)
                type = 'Vec';
                outstr = sprintf('%s; d_10 = %.2g, m_3 = %.2g',...
                    type,O.moments(1)/O.moments(0),O.moments(3) );
            else
                outstr = 'Empty';
            end % if else
            
            
        end
        
        %% - plot
        
        function pl_handle = plot(O,Parent)
            
            % Distribution.plot
            %
            % Plots number- and volume-weighted distributions
            %
            % SEE ALSO
            % Distribution
            
            if nargin < 2 || all(isempty(Parent)) || ~any(ishandle(Parent))
                Parent = figure;
                set(Parent,'numbertitle','off','name','PSDs (overlapping)');
            end % if
            
            Fax(1) = subplot(1,2,1,'Parent',Parent);
            Fax(2) = subplot(1,2,2,'Parent',Parent);
            xlabel(Fax(1),'Mean Char. Length')
            xlabel(Fax(2),'Mean Char. Length')
            ylabel(Fax(1),'Normalized Number Distribution')
            ylabel(Fax(2),'Normalized Volume Distribution')
            
            box(Fax(1),'on');
            box(Fax(2),'on');
            
            hold(Fax(1),'all')
            hold(Fax(2),'all')
            
            pl_handle = zeros(2*length(O),1);
            
            for i = 1:length(O)
                pl_handle(2*i-1) = plot(O(i).y,O(i).F./moments(O(i),0),'Parent',Fax(1));
                pl_handle(2*i) = plot(O(i).y,O(i).F.*O(i).y.^3./moments(O(i),3),'Parent',Fax(2));
            end % for
            
            if nargout < 1
                clear pl_handle
            end % if
            
        end % function
        
        %% - isnan
        
        function no = isnan(O)
            
            % Distribution.isnan
            %
            % Always returns false
            %
            % SEE ALSO
            % Distribution
            
            no = false;
        end
        
    end % methods
    
end % classdef