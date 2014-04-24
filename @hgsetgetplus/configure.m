function config = configure(obj)
% CONFIGURE Return validation criteria for SET.
%
%   CONFIG = CONFIGURE(OBJ) inputs an HGSETGETPLUS object and
%   returns a structure with one field per property in OBJ.
%   Each field is itself a structure with all of the following
%   fields:
%
%       classes: cell array of strings
%       attributes: cell array, mixed string and numeric
%       choices: cell array of strings
%       default: string
%
%   The method CONFIGUREDEFAULT initializes CONFIG so all
%   properties are present with all of their fields and the
%   default values. If any are missing, use of the class may
%   return errors.
%
%   For allowable classes, attributes, and choices,
%   see also: VALIDATEATTRIBUTES, VALIDATESTRING

% Keep attributes stored to avoid repeated calculations.
persistent S

if isempty(S)
    
    % Keep this line! The structure needs to be initialized.
    S = configureDefault(obj);
    
    % Now insert any constraints that are desired in one of
    % four fields: classes, attributes, choices, default. All
    % fields are cell arrays except default, which is a string.
    %
    % Example:
    %  S.Prop1.classes = {'char'};
    %  S.Prop1.choices = {'SI', 'cgs'};
    %  S.Prop1.default = 'SI';
end

config = S;