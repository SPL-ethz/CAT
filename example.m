%% Clean up

clear
clc
close all

%% Set up problem

% Basic object
PD = ProblemDefinition

% Define initial conditions
PD.init_dist.y = linspace(0,2);
PD.init_dist.F = @(x) normalpdf(x,0.5,0.1,10);
PD.init_conc = 1;

% Define growth rate
PD.growthrate = @(c,y)1e-2*ones(size(y));

% Define solution parameters - time
PD.sol_time = [0 100];

%% Solve

[t_out SolF Solc] = PBESolver(PD);

%% Plot results

% Plot distributions
pls = plot(SolF);

% Plot concentration
figure
plot(t_out,Solc)
xlabel('Time')
ylabel('Concentration')