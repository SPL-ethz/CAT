%% Clean up

clear all
clc
close all

addALLthepaths


%% Set up problem

% Basic object
PD = ProblemDefinition;

% Define initial conditions
nBins = 50;
gridL = linspace(0,5e2,nBins+1);
meanL = (gridL(1:end-1)+gridL(2:end))/2;
PD.init_dist.y = meanL;
PD.init_dist.boundaries = gridL;

PD.init_dist.F = @(x) normpdf(x,100,30);
PD.init_seed = 1;
PD.init_massmedium = 1000;

PD.sol_method = 'movingpivot';

% Define growth rate
PD.growthrate = @(S,T,y) 1e-1*log(S)*ones(size(y));

% Define solution parameters - time
PD.sol_time = [0 60*60];

%% Solve

PD = ProfileManager(PD);


%% Plot results

plot(PD,'detailed_results');