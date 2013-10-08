%% Clean up

if ~exist('demo','var')
    addALLthepaths
    clear all
    clc
    close all
    
    PD = ProblemDefinition;
    
%     PD.sol_method = 'movingpivot';
%     PD.sol_method = 'centraldifference';
    PD.sol_method = 'hires';
    nBins = 100;
end



%% Set up problem

% Define grid
% 
gridL = linspace(0,5e2,nBins+1);
meanL = (gridL(1:end-1)+gridL(2:end))/2;
PD.init_dist.y = meanL;
PD.init_dist.boundaries = gridL;
PD.init_conc = 6.5e-3;

% Define growth rate
PD.growthrate = @(S,T,y) 1e-1*(S-1)*ones(size(y));
PD.nucleationrate = @(S,T,F) (S>1)*exp(-100/log(S)/T)*moments(F,3)/moments(PD.init_dist,3);
% PD.nucleationrate = @(S,T) (S>1)*exp(-100/log(S)/T);
PD.solubility = @(T) (0.0056*(T-273).^2+0.0436.*(T-273)+3.7646)/1000;

% Define operating conditions
PD.init_seed = 1;
PD.init_massmedium = 2000; % mass of solvent in the beginning
PD.sol_time = [0 60*60];
PD.Tprofile = [0 5*60 10*60 60*60;
    290 285 285 280];

%define a simple gaussian as initial distribution
mu = 1e2;
sigma = 0.3*mu;
gauss = @(x) exp(-((x-mu).^2/(2*sigma^2)));

PD.init_dist.F = gauss(meanL);


PD = ProfileManager(PD);



%% Plot results
if ~exist('demo','var')
    plot(PD,'detailed_results');
end