function configStr = configToString(config)
% Convert a structure containing validation information to a structure
% containing strings.
%
% CONFIGSTR = CONFIGTOSTRING(INSTRUCT) inputs INSTRUCT, a structure
% containing validation information for a property. Returns CONFIGSTR, a
% string representing the information.
%
% See HGSETGETPLUS.CONFIGURE for a description of the contents of the input
% structure. This method does not check for correctness of the input. If
% all of the attributes are empty, returns an empty string.
%
% There is a hierarchy of constraints. If CHOICES is not empty, then the
% only other field that is read is DEFAULT. If CHOICES is empty but
% ATTRIBUTES is nonempty, then CLASSES and DEFAULT are ignored. Only if
% CHOICES and ATTRIBUTES are empty then CLASSES is read.


configStr = '';

if ~isempty(config.choices)
    in = config.choices;
    
    % Find default, if any, and put brackets around it.
    if ~isempty(config.default) && ischar(config.default)
        
        indx = find(strcmp(config.default,config.choices));
        if any(indx)
            in{indx} = ['{', in{indx}, '}'];
        end
    end
    configStr = sprintf('| %s ', in{:});
    configStr = ['[', configStr(2:end), ']'];
    return
end
    
if ~isempty(config.classes)
    in = config.classes;
    
    configStr = in{1};
    for j=2:length(in)
        configStr = [configStr,' -or- ',in{j}];
    end
end

if ~isempty(config.attributes)
    in = config.attributes;
    
    if length(config.classes) > 1
        configStr = ['(',configStr,')'];
    end
    if ~isempty(configStr)
        configStr = [configStr,' -and- '];
    end
    
    configStr = [configStr,in{1}];
    for j=2:length(in)
        if isnumeric(in{j})
            configStr = [configStr,' ',mat2str(in{j})];
        else
            configStr = [configStr,' -and- ',in{j}];
        end
    end
end