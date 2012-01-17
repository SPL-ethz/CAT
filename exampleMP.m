%% Clean up

clear all
clc
close all

%% Set up problem

% Basic object
PD = ProblemDefinition;

% Define grid
nBins = 100;
gridL = linspace(0,1000e-6,nBins+1);
meanL = (gridL(1:end-1)+gridL(2:end))/2;
PD.init_dist.y = meanL;
PD.init_dist.boundaries = gridL;
PD.init_conc = 3*SolubilityAlphaLGLU(PD.init_temp);

% Define growth rate
PD.growthrate = @(c,T,y) GrowthRateAlphaLGLU(c/SolubilityAlphaLGLU(T),T,y);
% PD.growthrate = @(c,T,y) zeros(1,length(y));

% Define nucleation rate
PD.nucleationrate = @(c,T) NucleationRateAlphaLGLU(c/SolubilityAlphaLGLU(T),T); %comment to deactivate nucleation

% Define operating conditions
seed_mass = 0.004; % seed mass - kg
PD.init_volume = 0.02; % volume of reactor - m^3
PD.sol_time = linspace(0,10000,51); % time the process is run
% PD.sol_time = [0,3600]; % time the process is run
PD.coolingrate = -0.0069; % in [K/s], equals 25�C/hr

%define a simple gaussian as initial distribution
mu = mean(gridL);
sigma = 0.1*mu;
gauss = @(x) exp(-((x-mu).^2/(2*sigma^2)));
values = gauss(meanL);
values = values.*(gridL(2:end)-gridL(1:end-1)); %transforming it into a number distribution 
values = seed_mass/(PD.rhoc*PD.kv*sum(values.*meanL.^3))/PD.init_volume*values(:);
values = values./(gridL(2:end)-gridL(1:end-1))';
PD.init_dist.F = values;


% Set solver method to moving pivot
PD.sol_method = 'movingpivot';


%% Solve
[PD.calc_time, PD.calc_dist, PD.calc_conc, PD.calc_temp, PD.calc_volume] = PBESolver(PD);

%% Plot results

plot(PD,'detailed_results');
