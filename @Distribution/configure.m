function config = configure(obj)
% CONFIGURE Return validation criteria for SET.
%
%   CONFIGSTRUCT = CONFIGURE(OBJ) inputs an HGSETGETPLUS object
%   and returns a structure with one field per property in OBJ.
%   Each field is itself a structure with all of the following
%   fields:
%
%       classes: cell array of strings
%       attributes: cell array, mixed string and numeric
%       choices: cell array of strings
%       default: string
%
%   The method CONFIGUREDEFAULT initializes CONFIGSTRUCT so all
%   properties are present with all of their fields and the
%   default values. If any are missing, use of the class may
%   return errors.
%
%   For allowable classes, attributes, and choices,
%   see also: VALIDATEATTRIBUTES, VALIDATESTRING

% Defaults and possible choices are

% Keep attributes stored to avoid repeated calculations.
persistent DistributionConfig

if isempty(DistributionConfig)
    
    % Keep this line! The structure needs to be initialized.
    DistributionConfig = configureDefault(obj);
    
    % Now insert any constraints that are desired in one of
    % three fields: classes, attributes, and choices. The
    % default values should be assigned in the Properties
    % block.
    %
    
    DistributionConfig.y.classes = {'double'};
    DistributionConfig.y.attributes = {'row','nonempty','real','finite','nonnegative'};
    
    % Boundaries has same constraints as y
    DistributionConfig.boundaries = DistributionConfig.y;
    
    % Distribution: function or vector
    DistributionConfig.F.classes = {'double','function_handle'};
    DistributionConfig.F.attributes = {'nonnegative','real','finite','row'};
    
    % Mu and sigma
    DistributionConfig.mu.classes = {'double'};
    DistributionConfig.mu.attributes = {'scalar','real','finite'};
    
    DistributionConfig.sigma.classes = {'double'};
    DistributionConfig.sigma.attributes = {'scalar','real','finite','positive'};
    
    % Mass
    DistributionConfig.mass.classes = {'double'};
    DistributionConfig.mass.attributes = {'real','nonnegative','finite'};
       
end

config = DistributionConfig;

end % function