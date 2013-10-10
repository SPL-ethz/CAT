classdef Distribution < handle
    
    % Class Distribution
    %
    % Used for defining distributions.
    %
    % The size vector needs to be defined. The distribution can be defined
    % as either a function handle or as a vector. The class checks the
    % inputs, and returns a vector in any case.
    %
    % The class can be called as follows:
    %
    %   dist = Distribution(y,@(x)normalpdf(x,1,0.1));
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
    
    %% Properties
    
    properties
        
        % Initial pivot vector
        % Define as vector
        y = linspace(0,1,20);
        
        % Initial vector of boundaries
        % Define as vector
        boundaries
        
        
    end % properties
    
    properties (Dependent)
        
        % Initial distribution
        % Define as either a function or a vector
        F = [];
        
        % mass of distribution (in g)
        mass
        
    end % properties
    
    properties (Access=private)
        
        % Private version of F - this is where the data is stored.
        pF = [];
        
    end % properties
    
    methods
        
        %% Method Distribution (constructor)
        
        function O = Distribution(y,F,boundaries)
            
            % If given, define y and F based on inputs. Otherwise, do
            % nothing
            
            if nargin > 0 && ~isempty(y)
                O.y = y;
            end % if
            
            if nargin > 1 && ~isempty(F)
                O.F = F;
            end % if
            
            if nargin > 2 && ~isempty(boundaries)
                O.boundaries = boundaries;
            end % if
            
        end % function
        
        %% Method set.y
        
        function set.y(O,value)
            
            % SET.Y
            %
            % Check input for y: a vector, all positive, increasing. Always
            % set y as a row vector

            if all(value>=0) && all(isfinite(value)) && length(value) > 1 && isvector(value) && ~any(diff(value)<0)
                O.y = value(:)';
            else
                warning('Distribution:SetY:WrongValue',...
                    'The property y must be a vector of positive values, in increasing order (duplicates allowed)');
            end % if
            
        end % function
                    
        %% Method set.boundaries
        
        function set.boundaries(O,value)
            
            % SET.Y
            %
            % Check input for boundaries: a vector, all positive (and 0),
            % duplicates allowed. Always set boundaries as a row vector
            
            if all(value>=0) && all(isfinite(value)) && length(value) > 1 && isvector(value) && ~any(diff(value)<0)
                O.boundaries = value(:)';
            else
                warning('Distribution:SetBoundaries:WrongValue',...
                    'The property boundaries must be a vector of positive values, in increasing order (duplicates allowed)');
            end % if
            
        end % function
        
        %% Method set.F
        
        function set.F(O,value)
            
            % SET.F
            %
            % Check input: vector, function handle, empty
            
            if isa(value,'function_handle')
                % Value is ok, set                
                    O.pF = value;                
            elseif isvector(value) || isempty(value)
                O.pF = value(:)'; %make it a row vector
            elseif iscell(value) && length(value)==3
                if strcmpi(value{1},'normal')
                   O.pF = @(x) 1/(value{3}*sqrt(2*pi))*exp(-((x-value{2}).^2/(2*value{3}^2))); 
                elseif strcmpi(value{1},'lognormal')
                   O.pF = @(x) 1/(x*value{3}*sqrt(2*pi))*exp(-((log(x)-value{2}).^2/(2*value{3}^2))); 
                end
            else
                % Value is not OK - display warning
                warning('Distribution:setF0:WrongType',...
                    'F has to be a function handle, a vector or a cell of length 3');
            end % if
            
        end % function
        
        %% Method get.F
        
        function Fout = get.F(O)
            
            % GET.F
            % Return F, based on what is in y.
            % If F is a function handle: calculate vector using positions
            % at y. If F is a vector, check so that its size is the same
            % as that of y
            
            if isa(O.pF,'function_handle')
                % Function handle - calculate values
                
                Fout = O.pF(O.y);
                
            elseif isvector(O.pF)
                
                if any( size(O.pF) ~= size(O.y) )
                    warning('Distribution:getF0:WrongSize',...
                        'F is not the same size as y');
                    Fout = [];
                else
                    Fout = O.pF;
                end % if else
                
            else
                Fout = [];
            end % if else
            
            % Always return a row vector
            Fout = Fout(:)';
            
        end % function
        
        %% Method set.mass
        
        function set.mass(O,value)
            
            O.F = O.F*value(1)/(moments(O,3)*value(2)*value(3)*value(4));
                        
        end % function
        
        %% Method moments(F,j)
        
        function Fmo = moments(O,j,icalc)
            
            % Calculate moments from distribution/distributions over time.
            %
            % Use moments(F,j) to calculate the j'th moment of F over its
            % length. 
            %
            % moments returns the indicated moment as a row vector
            %
            % To calculate moments at specific time points, use e.g.:
            %   F.moments(3,[1 3 5])
            % This returns the moment at the times equal t_out([1 3 5])
            
       
            if nargin <= 2 && ~exist('icalc','var')
                icalc = 1:length(O);
            end
            
            % If possible, calculate moment. Otherwise do nothing
            if nargin > 1 && ~isempty(j) 
                
                % Preallocate Fmo
                Fmo = zeros(size(icalc));
                
                for i = icalc
                    if isempty(O(1).boundaries)
                        Dy = diff([0 O(i).y]);
                    else
                        Dy = diff(O(i).boundaries);
                    end
                    
                    Fmo(icalc==i) = sum(O(i).F(:) .* Dy(:).* O(i).y(:).^j);
                end % for
                    
            else
                warning('Distribution:moments:nomoment',...
                    'No type of moment was indicated');
                Fmo = [];
            end % end
            
            % Always return a row vector
            Fmo = Fmo(:)';
            
        end % function

        %% Method plot
        
        function pl_handle = plot(O)
            
            nargout
            
            % Plot distribution, number- and volume-weighted
            
            FFig = figure;
            set(FFig,'numbertitle','off','name','PSDs (overlapping)')
            Fax(1) = subplot(1,2,1);
            Fax(2) = subplot(1,2,2);
            xlabel(Fax(1),'Mean Char. Length')
            xlabel(Fax(2),'Mean Char. Length')
            ylabel(Fax(1),'Number Distribution')
            ylabel(Fax(2),'Normalized Volume Distribution')
            
            box(Fax(1),'on')
            box(Fax(2),'on')
            
            hold(Fax(1),'all')
            hold(Fax(2),'all')
            
            pl_handle = zeros(2*length(O),1);
            
            for i = 1:length(O)
                pl_handle(2*i-1) = plot(O(i).y,O(i).F,'Parent',Fax(1));
                pl_handle(2*i) = plot(O(i).y,O(i).F.*O(i).y.^3,'Parent',Fax(2));
            end % for
            
        end % function

    end % methods
    
end % classdef