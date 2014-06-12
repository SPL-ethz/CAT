function setDefaults(Oin,varargin)

% setDefaults
%
% set sensible defaults for CAT method
%

for O = Oin
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.init_dist)
        O.init_dist = Distribution(linspace(1,750,250),{'normal',100,20});
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.Tprofile)
        O.Tprofile = @(t) 25;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.ASprofile)
        O.ASprofile = @(t) 0*t;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.init_conc)
        O.init_conc = 1.2; % supersaturated
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.solubility)
        O.solubility = @(T) 1;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.growthrate)
        O.growthrate = @(S) (S-1);
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.nucleationrate)
        O.nucleationrate = @(S) 0;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.init_seed)
        O.init_seed = 1.5;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.init_massmedium)
        O.init_massmedium = 1000;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.rhoc)
        O.rhoc = 1e-12;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.kv)
        O.kv = 1;
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.sol_method)
        O.sol_method = 'cd';
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.sol_method)
        O.sol_options= {''};
    end
    
    if ~any(strcmp(varargin,'emptyonly')) || isempty(O.sol_time)
        O.sol_time = 0:100:10000;
    end
    
end % for

end