%% Clean up

clear all
clc
close all

addALLthepaths


%% Set up problem

% Basic object
PD = ProblemDefinition;

% Define initial conditions
PD.init_dist.y = linspace(10,500);
PD.init_dist.F = @(x) normpdf(x,100,30);
PD.init_seed = 1;
PD.init_massmedium = 1000;

% Define growth rate
PD.growthrate = @(S,T,y) 1e-1*log(S)*ones(size(y));

% Define solution parameters - time
PD.sol_time = [0 60*60];

%% Solve

PD = ProfileManager(PD);


%% Plot results

plot(PD,'detailed_results');