%% Clean up

if ~exist('demo','var')
    clear all
    clc
    close all
    
    cat = CAT;
    
%     cat.sol_method = 'movingpivot';
%     cat.sol_method = 'centraldifference';
    cat.sol_method = 'hires';
    nBins = 100;
end

% Define grid
gridL = linspace(0,5e2,nBins+1);
meanL = (gridL(1:end-1)+gridL(2:end))/2;
cat.init_dist.y = meanL;
cat.init_dist.boundaries = gridL;
cat.init_conc = 7e-3;

% Define growth rate
cat.growthrate = @(S,T,y) (S>1)*4e-1*(S-1)*ones(size(y));
cat.solubility = @(T) (0.0056*(T-273).^2+0.0436.*(T-273)+3.7646)/1000;

% Define operating conditions
cat.init_seed = 0.5;
cat.init_massmedium = 2000; % mass of solvent in the beginning
cat.sol_time = [0 60*60];
cat.Tprofile = [0 5*60 10*60 60*60;
    290 285 285 280];
cat.ASprofile = [0 5*60 10*60 60*60;
    0 50 200 200];

%define a simple gaussian as initial distribution
mu = 1e2;
sigma = 0.3*mu;
gauss = @(x) exp(-((x-mu).^2/(2*sigma^2)));

cat.init_dist.F = gauss(meanL);


% Set solver method to moving pivot
% cat.sol_method = 'movingpivot';
% cat.sol_method = 'centraldifference';
% cat.sol_method = 'hires';
cat = ProfileManager(cat);



%% Plot results
if ~exist('demo','var')
    plot(cat,'detailed_results');
end