classdef Easyset < hgsetget
    %% Easyset
    % Class to implement some automation on set methods
    %
    % Check set value against defined:
    % - choices of variables
    % or
    % - list of classes and list of attributes
    %
    % How to use this class:
    %
    % Make your class a subclass of Easyset:
    %   classdef NewClass < Easyset
    %
    % Define set methods for each variable which you want to restrict. The
    % class should:
    %  1. Define the choices/classes/attributes structure for the current
    %  property
    %  2. Call the setPropertyValue method to perform the checking (will
    %  call error if any checks fail)
    %  3. Do any optional extra setting
    %  4. Save the value at the end of the function
    %
    %% EXAMPLE for a property named VAR:
    %
    %     function set.VAR(O,value)
    %
    %         % Define properties here
    %         O.classes.VAR = 'numeric';
    %         O.attributes.VAR = {'vector','real','finite','nonnegative'};
    %
    %         % Redirect to setPropertyValue function to do checking
    %         O.setPropertyValue('VAR',value);
    %
    %         % If checks didn't fail, method has not exited, so do extra
    %         % checks here
    %
    %         % Set the value
    %         O.VAR = value;
    %
    %     end % function
    %
    % Instead of defining the property structure at the beginning of the
    % set method, the method configureproperties can be defined to take
    % care of this. This can be useful when many properties share the same
    % definition, for example VAR and VAR2 in this example:
    %
    %     function configureproperties(O)
    %
    %         % In this case, VAR and boundaries have the same restrictions
    %         [O.classes.VAR,O.classes.VAR2] = deal({'numeric'});
    %         [O.attributes.VAR,O.attributes.VAR2] = deal({'vector','real','finite','nonnegative'});
    %
    %     end % function
    %
    %% SEE ALSO
    % Easyset.choices, Easyset.classes, Easyset.attributes,
    % Easyset.configureproperties

    %% Easyset properties
    
    properties (Access=protected)
        
        % Easyset.choices
        %
        % Structure of choices for restricted properties. These override
        % any attribute and class choices. Defined as a structure listing
        % the allowable choices for each restricted property
        %
        % SEE ALSO
        % Easyset, Easyset.classes, Easyset.attributes
        choices
        
        % Easyset.classes
        %
        % Structure of allowed classes for restricted properties. Defined
        % as a structure with fields for each property, containing cell
        % values listing allowable classes.
        %
        % At least one class must be defined for each property if
        % attributes are to be checked! If no class is defined, no checking
        % is done on the attributes of the property.
        %
        % SEE ALSO
        % Easyset, Easyset.choices, Easyset.attributes
        classes
        
        % Easyset.attributes
        %
        % Structure of attributes for restricted properties. Defined as a
        % structure with fields for each property, containing cell values
        % listing allowable attributes.
        %
        % Only works if classes is defined for the property!
        %
        % SEE ALSO
        % Easyset, Easyset.choices, Easyset.classes
        attributes
        
    end
    
    %% Easyset methods
    methods
        
        %% Constructor method
        
        function O = Easyset()
            
            % Call configureproperties method - to be overwritten
            O.configureproperties;
            
        end % function
        
        %% Method configureproperties
        
        function configureproperties(O)
            
            % configureproperties
            %
            % This method is run when Easyset is constructed. It can be
            % used to define the structures choices, classes and attributes
            %
            % SEE ALSO
            % Easyset
            
            % Empty - to be overwritten by subfunctions
            
        end % function
        
        %% General set method
        % Call another method, depending on the number of inputs to this
        % function. Checking and settings is done by set.property functions
        % below
        
        function set(O,property,value)
            
            % set()
            % Set the value of a property:
            %   set(obj,property,value)
            % or get the values which can be set for that property:
            %   set(obj,property)
            %
            % SEE ALSO
            % Easyset
            
            % nargin: remember that O counts as well
            
            if nargin == 3 && ~isempty(property)
                % Attempt to set the value - use its set method
                O.(property) = value;
            elseif nargin == 2 && ~isempty(property)
                % Set without any values - send this to setPropertyValue so
                % that it can do its thing. Don't set anything
                O.setPropertyValue(property);
            else
                % Don't know what to do, show help for this function
                help('Distribution.set')
            end % if elseif else
            
        end % function
        
        %% Function setPropertyValue
        
        function setPropertyValue(O,property,varargin)
            
            % setPropertyValue
            % Function to check values against structure of choices,
            % classes and attributes defined in class file, or to return a
            % string describing what a property can be set to.
            %
            % SEE ALSO
            % Easyset
            
            % Decide what to do
            if length(varargin) == 1
                % Value requested to be set. Check it against the allowed
                % choices, classes and attributes
                
                value = varargin{1};
                
                % Make sure the property exists
                if isprop(O,property)
                    
                    % Check if limited choices are set - if not, then
                    % classes and attributes are checked
                    if isfield(O.choices,property) && ~isempty(O.choices.(property))
                        if ~any(strcmp(value,O.choices.(property)))
                            error(['Easyset:set:' property],...
                                'The parameter %s must be one of the following: %s.',...
                                property,implode(O.choices.(property),', '));
                        end % if
                        
                    else
                        
                        % No list of choices - check classes and attributes
                        
                        % Check if any class limitations are set
                        if isfield(O.classes,property) && ~isempty(O.classes.(property))
                            
                            if ~iscell(O.classes.(property))
                                O.classes.(property) = {O.classes.(property)};
                            end % if
                            
                            % Check if attribute limitations are set
                            if isfield(O.attributes,property) && ~isempty(O.attributes.(property))
                                checkattributes = O.attributes.(property);
                            else
                                checkattributes = {''}; % No attributes need to be defined
                            end % if
                            
                            % Do the check
                            validateattributes(value,O.classes.(property),checkattributes,'',property);
                            
                        end % if
                        
                    end % if
                    
                else
                    error('Easyset:set:invalidproperty',...
                        '%s is not a valid property',property);
                end % if else
                
            elseif isempty(varargin)
                % No values given to set to - return a string showing what
                % the property classes and attributes are
                
                % Initialise the string
                outstr = sprintf('Property %s ',property);
                
                % Make sure the property exists
                if isprop(O,property)
                    
                    % Check if limited choices are set - if not, then
                    % classes and attributes are checked
                    if isfield(O.choices,property) && ~isempty(O.choices.(property))
                        outstr = sprintf('%s has be one of the following values: %s.',outstr,implode(O.choices.(property),', '));
                    else
                        
                        % Check if any class limitations are set
                        if isfield(O.classes,property) && ~isempty(O.classes.(property))
                            outstr = sprintf('%s has to be of class %s;',...
                                outstr,implode(O.classes.(property),', '));
                        else
                            outstr = sprintf('%s has no class limitations;',outstr);
                        end % if
                        
                        % Check if attribute limitations are set
                        if isfield(O.attributes,property) && ~isempty(O.attributes.(property))
                            outstr = sprintf('%s has to have the following attributes: %s.',...
                                outstr,implode(O.attributes.(property),', '));
                        else
                            outstr = sprintf('%s has no attribute limitations.',outstr);
                        end % if
                        
                    end % if
                    
                else
                    outstr = sprintf('%s is not a valid property',property);
                end % if else
                
                % Display the output string
                fprintf('%s\n',outstr);
                
            else
                error('Easyset:set:toomanyinputs','Too many inputs for set');
            end % if elseif else
            
        end % function
        
    end
    
end