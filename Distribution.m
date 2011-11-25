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
    %   dist.F = @(x) normalpdf(x,1,0.1);
    %
    % 
    
    %% Properties
    
    properties
        
        % Initial grid vector
        % Define as vector
        y = linspace(0,1,20);
        
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
        
        function O = Distribution(y,F)
            
            % If given, define y and F based on inputs. Otherwise, do
            % nothing
            
            if nargin > 0 && ~isempty(y)
                O.y = y;
            end % if
            
            if nargin > 1 && ~isempty(F)
                O.F = F;
            end % if
            
        end % function
        
        %% Method set.F
        
        function set.F(O,value)
            
            % SET.F
            %
            % Check input: vector, function handle, or empty
            
            if strcmp(class(value),'function_handle') || isvector(value) || isempty(value)
                % Value is ok, set
                
                O.pF = value;
                
            else
                % Value is not OK - display warning
                warning('Distribution:setF0:WrongType',...
                    'F has to be a function handle or a vector');
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
            
        end % function
        
        %% Method plot(F)
        
        function plot(O,varargin)
            
            plot(O.y,O.F,varargin{:})
            xlabel('Size y')
            ylabel('Distribution F')
            
        end % function
        
    end % methods
    
end % classdef