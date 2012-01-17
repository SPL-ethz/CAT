%% Clean up

clear all
clc
close all

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
tprofile = [0 3600 21600 84000];
Tprofile = [308 303 298 298];


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
PD.sol_method = 'movingpivot';


%% Solve
PD.calc_dist = Distribution;
for i = 2:length(Tprofile)
    i
    PD.coolingrate = (Tprofile(i)-Tprofile(i-1))/(tprofile(i)-tprofile(i-1));
    PD.sol_time = [tprofile(i-1) tprofile(i)];  
    [a b c d e] = PBESolver(PD);
    PD.calc_time(end+1:end+length(a)) = a;
    PD.calc_dist(end+1:end+length(b)) = b;    
    PD.calc_conc(end+1:end+length(a)) = c; 
    PD.calc_temp(end+1:end+length(a)) = d;
    PD.calc_volume(end+1:end+length(a)) = e;
 
    PD.init_dist = PD.calc_dist(end);
    PD.init_temp = PD.calc_temp(end);
    PD.init_volume = PD.calc_volume(end);
    PD.init_conc = PD.calc_conc(end);
end

PD.calc_dist = PD.calc_dist(2:end);
PD.init_temp = Tprofile(1);
PD.sol_time = [tprofile(1) tprofile(end)];
PD.init_dist = PD.calc_dist(1);
% PD.seed_mass = seed_mass;



%% Plot results
plot(PD,'detailed_results');