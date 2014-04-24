classdef hgsetgetplus < hgsetget
    % HGSETGETPLUS   HG-style set and get for MATLAB objects.
    %     The hgsetgetplus class is an abstract class that provides an HG-style
    %     property set and get interface.  hgsetgetplus is a subclass of handle, so
    %     any classes derived from hgsetgetplus are handle classes.
    %
    %     classdef MyClass < hgsetgetplus makes MyClass a subclass of hgsetgetplus.
    %
    %     Classes that are derived from hgsetgetplus inherit no properties but
    %     do inherit methods that can be overridden as needed.
    %
    % HGSETGETPLUS Methods:
    %   set - set MATLAB object property values.
    %   get - get MATLAB object property values.
    %   setdisp - Specialized MATLAB object property display.
    %   getdisp - Specialized MATLAB object property display.
    %   configure - Return validation criteria for SET.
    %
    % See also HANDLE, HGSETGET
    
    methods
        function varargout = set(obj, varargin)
            % SET set MATLAB object property values.
            %     SET(H,'PropertyName',PropertyValue) sets the value of the
            %     specified property for the MATLAB object with handle H.
            %     If H is an array of handles, the specified property's
            %     value is set for all objects in H.
            %
            %     SET(H,'PropertyName1',Value1,'PropertyName2',Value2,...)
            %     sets multiple property values with a single statement.
            %
            %     SET(H,pn,pv) sets the named properties specified in the
            %     cell array of strings pn to the corresponding values in
            %     the cell array pv for all objects specified in H.  The
            %     cell array pn must be a vector of length N, but the cell
            %     array pv can be 1-by-N or M-by-N, where M is equal to
            %     length(H). If pv is 1-by-N, each property in pn is set to
            %     the same value in all objects in H. If pv is M-by-N, each
            %     object will be updated with a different set of values for
            %     the list of property names contained in pn.
            %
            %     Given S, a structure array whose field names are object
            %     property names, SET(H,S) sets the properties identified
            %     by each field name of S with the values contained in the
            %     structure. S can have length 1 or M, where M is equal to
            %     length(H). If S has length 1, each field name in S sets
            %     the corresponding property for all objects in H. If S has
            %     length M, each object will be updated with a different
            %     set of values for the list of property names contained in
            %     pn.
            %
            %     A = SET(H, 'PropertyName') returns the possible values
            %     for the specified property of the object with handle H.
            %     The returned array is a cell array of possible value
            %     strings or an empty cell array if the property does not
            %     have a finite set of possible string values.
            %
            %     SET(H,'PropertyName') displays the possible values for
            %     the specified property of object with handle H.
            %
            %     Note that it is permissible to use property/value string
            %     pairs, structures, and property/value cell array pairs in
            %     the same call to SET.
            %
            %     A = SET(H) returns all property names and their possible
            %     values for the object with handle H.  H must be scalar.
            %     The return value is a structure whose field names are the
            %     property names of H, and whose values are cell arrays of
            %     possible property values or empty cell arrays.
            %
            %     SET(H) displays all properties and property values of
            %     scalar object H.
            %
            %     The class can override the method SETDISP to control how
            %     information is displayed in SET(H), A = SET(H),
            %     SET(H,'PropertyName') and A = SET(H,'PropertyName').
            %
            %     See also: HGSETGETPLUS/CONFIGURE.
            
            
            if nargin == 1 % set(h) syntax
                
                
                if length(obj)>1
                    error('MATLAB:class:ScalarObjectRequired', ...
                        'Object passed to method ''set'' must be scalar.')
                end
                if nargout
                    varargout{1} = structfun(@(x) getfield(x,'choices'), ...
                        configure(obj),'UniformOutput',false);
                else
                    setdisp(obj);
                end
                
            elseif nargin == 2 && ~isstruct(varargin{1}) % set(h, 'PropertyName') syntax
                
                
                if length(obj)>1
                    error('MATLAB:class:ScalarObjectRequired', ...
                        'Object passed to method ''set'' must be scalar.')
                end
                np = varargin{1};
                
                if isempty(np) || (~iscell(np) && ~ischar(np))
                    error('MATLAB:class:InvalidArgument', ...
                        ['Invalid input argument type to ''set''.', ...
                        ' Type ''help hgsetgetplus/set'' for options.'])
                    
                elseif ~ismember(np, properties(obj))
                    error('MATLAB:class:PropertyNotFound', ...
                        ['Property %s not found in class %s,',...
                        ' or is not present in all elements of the array', ...
                        ' of class %s.'], ...
                        np, class(obj),class(obj));
                end
                
                if nargout
                    configStruct = configure(obj);
                    varargout{1} = configStruct.(np).choices;
                    
                else
                    setdisp(obj,np);
                end
                
            else % loop through variables
                
                iarg = 1;
                
                while iarg < nargin
                    p = varargin{iarg};
                    
                    if isstruct(p) % set(h,S) syntax: convert to cell
                        pn = fieldnames(p);
                        if numel(pn) == 1
                            pv = squeeze(struct2cell(p));
                        elseif ~isempty(pn)
                            pv = squeeze(struct2cell(p))';
                        end
                    else
                        iarg = iarg+1;
                        if iarg >= nargin % ran out of arguments
                            error('MATLAB:class:BadParamValuePairs', ...
                                'Invalid parameter/value pair arguments')
                        end
                        pn = p;
                        pv = varargin{iarg};
                    end
                    
                    if iscell(pn) % set(h, pn, pv) syntax
                        
                        if isempty(pn)
                            error('MATLAB:class:InvalidArgument', ...
                                ['Invalid input argument type to ''set''.', ...
                                ' Type ''help set'' for options.'])
                        end
                        if size(pn,1) > 1 && size(pn,2) > 1
                            error('MATLAB:class:ParamValueCellMismatch', ...
                                'Parameter must be scalar or vector.')
                        end
                        pn = pn(:)';
                        
                        
                        if length(obj)>1 % deal with arrays recursively
                            
                            for i=1:length(obj)
                                if size(pv,1) <= 1
                                    set(obj(i),pn,pv);
                                elseif size(pv,1) == length(obj)
                                    set(obj(i),pn,pv(i,:));
                                else
                                    error('MATLAB:class:ValueCellDimension', ...
                                        ['Value cell array handle dimension',...
                                        ' must match handle vector length.'])
                                end
                            end
                            
                        elseif size(pn,2) ~= size(pv,2)
                            error('MATLAB:class:ParamValueCellMismatch', ...
                                'Size mismatch in Param Cell / Value Cell pair.')
                            
                        else
                            for i=1:length(pn)
                                setOneProperty(obj,pn{i},pv{i});
                            end
                        end
                        
                        
                    elseif ischar(pn) % set(h, 'PropertyName', value) syntax
                        
                        for i=1:length(obj)
                            setOneProperty(obj(i),pn,pv);
                        end
                        
                    else
                        error('MATLAB:class:BadParamValuePairs', ...
                            'Invalid parameter/value pair arguments.')
                    end
                    
                    iarg = iarg+1;
                end
                
            end
        end
        
        function setdisp(obj,np)
            % SETDISP Specialized MATLAB object property display.
            %   SETDISP is called by SET when SET is called in one of the
            %   following formats:
            %
            %   SET(H)
            %   SET(H,'PropertyName')
            %
            %   The arguments from SET are passed to SETDISP and the
            %   property information is displayed to the workspace.
            %
            %   This class is designed to be overridden if a specialized
            %   format is needed for displaying the class. If not
            %   overridden, SETDISP displays the class in the default
            %   format.
            %
            %   See also HGSETGETPLUS, HGSETGETPLUS/SET, HANDLE
            
            %   SETDISP uses the hidden static method CONFIGTOSTRING to
            %   produce the string for each property.
            
            config = configure(obj);
            
            if nargin == 1
                
                names = fieldnames(config);
                outStruct = cell2struct(cell(size(names)),names);
                
                for i=1:length(names)
                    np = names{i};
                    configString = obj.configToString(config.(np));
                    if isempty(configString)
                        outStruct.(np) = {};
                    else
                        outStruct.(np) = configString;
                    end
                end
                disp(outStruct)
                
            else
                configString = obj.configToString(config.(np));
                
                if isempty(configString)
                    configString = '     {}';
                end
                disp(configString)
            end
        end
    end
    
    methods(Access=protected,Hidden)
        setOnePropergy(obj,pname,pval)
        
        S = configureDefault(obj)
    end
    
    methods(Static, Hidden)
        configStr = configToString(config)
    end
end