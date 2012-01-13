%% Clean up

clear all
clc
close all

%% Set up problem

% Basic object
PD = ProblemDefinition;

% Define grid
nBins = 200;
gridL = linspace(0,3e-4,nBins+1);
meanL = (gridL(1:end-1)+gridL(2:end))/2;
PD.init_dist.y = meanL;
PD.init_dist.boundaries = gridL;
PD.init_conc = 5*SolubilityAlphaLGLU(PD.init_temp);

% Define growth rate
PD.growthrate = @(c,T,y) GrowthRateAlphaLGLU(c/SolubilityAlphaLGLU(T),T,y);

% Define nucleation rate
% PD.nucleationrate = @(c,T) NucleationRateAlphaLGLU(c/SolubilityAlphaLGLU(T),T); %comment to deactivate nucleation

% Define operating conditions
seed_mass = 0.004; % seed mass - kg
PD.init_volume = 0.02; % volume of reactor - m^3
PD.sol_time = linspace(0,1000,51); % time the process is run

%define a simple gaussian as initial distribution
mu = 5e-5;
sigma = 0.1*100e-6;
gauss = @(x) exp(-((x-mu).^2/(2*sigma^2)));
values = gauss(meanL);
values = seed_mass/(PD.rhoc*PD.kv*sum(values.*meanL.^3))/PD.init_volume*values(:);
PD.init_dist.F = values*1e6;


% Set solver method to moving pivot
PD.sol_method = 'hires';
PD.sol_options = {'Phi' 'vanleer'};

%% Solve
[PD.calc_time PD.calc_dist PD.calc_conc] = PBESolver(PD);

%% Plot results

plot(PD,'results');

rel_error=massbal(PD);