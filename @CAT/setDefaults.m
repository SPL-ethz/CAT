function setDefaults(Oin,varargin)

% setDefaults
%
% set sensible defaults for CAT method
%

% Copyright 2015-2016 David Ochsenbein
% Copyright 2012-2014 David Ochsenbein, Martin Iggland
% 
% This file is part of CAT.
% 
% CAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% 
% CAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


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