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

% Solubility function load from database
% Requires a working copy of mym 1.36
% Check http://129.132.152.27/PBEToolbox for administration and infos

solData.serv = '129.132.152.27';
solData.user = 'PBEToolbox';
solData.pass = 'toolbox2000';
solData.solSet = 1;

mym('open',solData.serv,solData.user,solData.pass);
mym('use',solData.user);
tmpString = mym('SELECT function FROM solubilityFunctions WHERE id = "{S}"',solData.solSet);
mym('close');

PD.solubility = eval(['@(T,t,x,L) ' char(cell2mat(tmpString.function))' ';']);

% Define growth rate
PD.growthrate = @(c,T,y) 1e-2*ones(size(y));

% Define solution parameters - time
PD.sol_time = [0 100];

%% Solve

[PD.calc_time PD.calc_dist PD.calc_conc] = PBESolver(PD);


%% Plot results

plot(PD,'results');