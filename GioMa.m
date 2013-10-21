
% This is the prototype template that is filled out by fillOutForm

kitty = CAT; 
t0 = 1000;
%--------------------------------------------------------------------------
%% SIMULATION SETTINGS
% Solution Time
% FIELD:    sol_time
% STATUS:   REQUIRED
% UNITS:    [s]
% CLASS:    VECTOR
% INPUTS:   -
% EXP.:     [0 60*60], [0:5*60:60]
% COMMENT:  When length(sol_time)>2, the result will only be given at the
%           times indicated in sol_time.
kitty.sol_time = [0 3600];

% Solution Method
% FIELD:    sol_method
% STATUS:   OPTIONAL
% UNITS:    -
% CLASS:    STRING
% INPUTS:   -
% EXP.:     'hires', 'movingpivot', 'centraldifference'
% COMMENT:  Default solution method is centraldifference.
kitty.sol_method = 'cd';

%--------------------------------------------------------------------------
%% THERMODYNAMICS & KINETICS
% Solubility 
% FIELD:    solubility
% STATUS:   REQUIRED
% UNITS:    [g/g]
% CLASS:    FUNCTION_HANDLE
% INPUTS:   1. Temperature [°C],  2.(optional) Mass fraction antisolvent [-])
% EXP.:     @(T,xm) 1e-4*(1-xm)+(0.0056*(T).^2+0.0436.*(T)+3.7646)/1000
% COMMENT:  -
kitty.solubility = @(T)(0.0056*(T).^2+0.0436.*(T)+3.7646)/1000;

% Growth Rate 
% FIELD:    growthrate
% STATUS:   REQUIRED
% UNITS:    [micron/s]
% CLASS:    FUNCTION_HANDLE
% INPUTS:   1. Supersaturation [-], 2. (optional) Temperature [°C], 3. (optional) Particle Length [micron]
% EXP.:     @(S,T,y) (S>1)*4e-1*(S-1)*ones(size(y))
% COMMENT:  -
kitty.growthrate = @(S,T,y)1e-1*(S-1)*ones(size(y));

% Nucleation
% Growth Rate 
% FIELD:    nucleationrate
% STATUS:   OPTIONAL
% UNITS:    [#/(g s)]
% CLASS:    FUNCTION_HANDLE
% INPUTS:   1. Supersaturation [-], 2. Temperature [°C], 3. (optional) Distribution
% EXP.:     @(S,T,F) exp(-100/log(S)/T)*moments(F,3)/moments(PD.init_dist,3)
% COMMENT:  -
kitty.nucleationrate = @(S,T,t) (t>=t0)*(S>1)*exp(-100/log(S)/(T+273));

%--------------------------------------------------------------------------
%% OPERATING CONDITIONS
% Seed mass
% FIELD:    init_seed
% STATUS:   OPTIONAL
% UNITS:    [g]
% CLASS:    SCALAR
% INPUTS:   -
% EXP.:     1
% COMMENT:  -
kitty.init_seed = 0;

% Initial Concentration 
% FIELD:    init_conc
% STATUS:   REQUIRED
% UNITS:    [g/g]
% CLASS:    SCALAR, STRING
% INPUTS:   -
% EXP.:     7e-3, 'sat'
% COMMENT:  -
kitty.init_conc = 'sat';

% Initial Seed Distribution
% FIELD:    init_dist
% STATUS:   REQUIRED
% UNITS:    []
% CLASS:    DISTRIBUTION (see Class Distribution on how to construct distributions)
% INPUTS:   -
% EXP.:     Distribution(linspace(1,500),{'normal',50,20})
% COMMENT:  -
kitty.init_dist = Distribution(linspace(1,500,100),... 
@(x)zeros(size(x)));

% Total initial mass of Solvent + Antisolvent
% FIELD:    massmedium
% STATUS:   REQUIRED
% UNITS:    [g]
% CLASS:    SCALAR
% INPUTS:   -
% EXP.:     2000
% COMMENT:  -
kitty.init_massmedium = 2000; 

% Temperature Profile 
% FIELD:    Tprofile
% STATUS:   OPTIONAL
% UNITS:    [°C]
% CLASS:    SCALAR, MATRIX, FUNCTION_HANDLE
% INPUTS:   1. Time [s]
% EXP.:     25, [0 5*60 10*60 60*60; 25 25 23 20], @(t) (33-31)/3600*t+33
% COMMENT:  When using function handles with non-smooth functions make sure
%           that nodes are in the solution time vector (sol_time)!
kitty.Tprofile = [0 300 600 3600;17 12 12 7];

% Antisolvent Mass Profile 
% FIELD:    ASprofile
% STATUS:   OPTIONAL
% UNITS:    [g]
% CLASS:    MATRIX, FUNCTION_HANDLE
% INPUTS:   1. Time [s]
% EXP.:     [0 5*60 10*60 60*60; 25 25 23 20], @(t) (33-31)/3600*t+33
% COMMENT:  The Mass Profile should be monotonically increasing!
%           When using function handles with non-smooth functions make sure
%           that nodes are in the solution time vector (sol_time)! 
kitty.ASprofile = [];



% SOLVE
kitty.solve;

%% PLOT (optional)
kitty.plot;
 
 % This file was generated for you by the fillOutForm function.