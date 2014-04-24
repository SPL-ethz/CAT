function setOneProperty(obj,pname,pval)
% Set a single property using a name/value pair.

assert(numel(obj)==1)

allprops = properties(obj);
if ~ischar(pname)
    error('MATLAB:class:InvalidArgument', ...
        ['Invalid input argument type to ''set''.', ...
        ' Type ''help hgsetgetplus/set'' for options.'])
end

if ~ismember(pname, allprops)
    error('MATLAB:class:InvalidProperty', ...
        ['The name "%s" is not an accessible property', ...
        ' for an instance of class "%s".'], ...
        pname, class(obj));
end

% Check validity of value
S = configure(obj);
Sfield = S.(pname);

% Choices must be a cell array of strings. If it is nonempty,
% make sure that the class is set to 'char'.
if ~isempty(Sfield.choices)
    Sfield.classes = {'char'};
end

% Check classes first
if ~isempty(Sfield.classes)
    validateattributes(pval,Sfield.classes,Sfield.attributes, ...
        mfilename,pname);
end

% Then check choices, if necessary.
if ~isempty(Sfield.choices)
    pval = validatestring(pval,Sfield.choices, ...
        mfilename, pval);
end

obj.(pname) = pval;