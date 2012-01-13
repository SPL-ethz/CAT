%% Clean up

clear all
clc
close all

%% Set up problem

% Basic object
PD = ProblemDefinition;

% Define initial conditions
PD.init_dist.y = linspace(0,2);
PD.init_dist.F = @(x) normpdf(x,0.5,0.1);

% Define growth rate
PD.growthrate = @(c,T,y) 1e-2*ones(size(y));

% Define solution parameters - time
PD.sol_time = [0 100];

%% Solve

[PD.calc_time, PD.calc_dist, PD.calc_conc, PD.calc_temp, PD.calc_volume] = PBESolver(PD);


%% Plot results

plot(PD,'detailed_results');