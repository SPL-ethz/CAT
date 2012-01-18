%% Clean up

clear all
clc
close all

addALLthepaths

%% Set up problem

% Basic object
PD = ProblemDefinition;

% Define grid
nBins = 100;
gridL = linspace(0,2.5e-4,nBins+1);
meanL = (gridL(1:end-1)+gridL(2:end))/2;
PD.init_dist.y = meanL;
PD.init_dist.boundaries = gridL;
PD.init_conc = 1*SolubilityAlphaLGLU(PD.init_temp);

% Define growth rate
PD.growthrate = @(c,T,y) GrowthRateAlphaLGLU(c/SolubilityAlphaLGLU(T),T,y);

% Define nucleation rate
% PD.nucleationrate = @(c,T) NucleationRateAlphaLGLU(c/SolubilityAlphaLGLU(T),T); %comment to deactivate nucleation

% Define operating conditions
seed_mass = 0.004;
PD.init_volume = 0.02; % volume of reactor - m^3
PD.tprofile = [0 3600 21600 84000];
PD.Tprofile = [290 295 295 290];
% PD.coolingrateprofile={@(t) -1e-6*t 0 1e-5};


%define a simple gaussian as initial distribution
mu = 5e-5;
sigma = 0.1*mu;
gauss = @(x) exp(-((x-mu).^2/(2*sigma^2)));
values = gauss(meanL);
values = values.*(gridL(2:end)-gridL(1:end-1)); %transforming it into a number distribution 
values = seed_mass/(PD.rhoc*PD.kv*sum(values.*meanL.^3))/PD.init_volume*values(:);
values = values./(gridL(2:end)-gridL(1:end-1))';
PD.init_dist.F = values;


% Set solver method to moving pivot
PD.sol_method = 'hires';


PD = ProfileManager(PD);



%% Plot results
plot(PD,'detailed_results');