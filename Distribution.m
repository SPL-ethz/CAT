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
            % Check input for y: a vector, all positive, monotonically
            % increasing. Always set y as a row vector
            
            if all(value>=0) && all(isfinite(value)) && length(value) > 1 && isvector(value) && ~any(diff(value)<0)
                O.y = value(:)';
            else
                warning('Distribution:SetY:WrongValue',...
                    'The property y must be a vector of positive values, in increasing order (duplicates allowed)');
            end % if
            
        end % function
        
        %% Method set.F
        
        function set.F(O,value)
            
            % SET.F
            %
            % Check input: vector, function handle, or empty
            
            if strcmp(class(value),'function_handle')
                % Value is ok, set                
                    O.pF = value;                
            elseif isvector(value) || isempty(value)
                O.pF = value(:)'; %make it a row vector
            else
                % Value is not OK - display warning
                warning('Distribution:setF0:WrongType',...
                    'F has to be a function handle or a vector');
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
        
        %% Method get.F
        
        function Fout = get.F(O)
            
            % GET.F
            % Return F, based on what is in y.
            % If F is a function handle: calculate vector using positions
            % at y. If F is a vector, check so that its size is the same
            % as that of y
            
            if strcmp(class(O.pF),'function_handle')
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
        
        %% Method moments(F,j)
        
        function Fmo = moments(O,j,icalc,varargin)
            
            % Calculate moments from distribution/distributions over time.
            %
            % Use moments(F,j) to calculate the j'th moment of F over its
            % length. 
            %
            % moments returns the indicated moment as a row vector
            %
            % To calculate moments at specific time points, use e.g.:
            % moments(F,3,[1 3 5])
            % This returns the moment at the times equal t_out([1 3 5])
            
       
            if nargin <= 2 && ~exist('icalc')
                icalc = [1:length(O)];
            end
            
            % If possible, calculate moment. Otherwise do nothing
            if nargin > 1 && ~isempty(j) 
                for i = icalc
                    Fmo(icalc==i) = sum(O(i).F .* O(i).y.^j);
                end % for
            else
                warning('Distribution:moments:nomoment',...
                    'No type of moment was indicated');
                Fmo = [];
            end % end
            
            % Always return a row vector
            Fmo = Fmo(:)';
            
        end % function

    end % methods
    
end % classdef