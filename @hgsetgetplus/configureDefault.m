function S = configureDefault(obj)
% Make sure that the attributes structure has all required
% fields.
S2 = struct;
S2.classes = {};
S2.attributes = {};
S2.choices = {};
S2.default = '';
names = properties(obj);
for i=1:length(names)
    S.(names{i}) = S2;
end

% If default values are assigned to the properties, make them
% the default
h = eval(class(obj));
for i=1:length(names)
    pval = get(h,names{i});
    if ~isempty(pval)
        S.(names{i}).default = pval;
    end
end
delete(h)