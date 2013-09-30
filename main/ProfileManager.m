function [PD] = ProfileManager(PD)

%% Solve
PD.calc_dist = Distribution;
% keyboard
PD.init_dist.mass = [PD.init_seed PD.kv PD.rhoc];

if ~isempty(PD.tNodes)
    for i = 2:length(PD.tNodes) % make sure you hit the different nodes of the non-smooth profiles
        
        PD.sol_time = [PD.tnodes(i-1) PD.tnodes(i)];  
        [a b c e] = PBESolver(PD);
        PD.calc_time(end+1:end+length(a)) = a;
        PD.calc_dist(end+1:end+length(b)) = b;    
        PD.calc_conc(end+1:end+length(a)) = c; 
        PD.calc_massmedium(end+1:end+length(a)) = e;

        PD.init_dist = PD.calc_dist(end);
        PD.init_conc = PD.calc_conc(end);

    end % for

        PD.calc_dist = PD.calc_dist(2:end);
        PD.init_temp = PD.calc_temp(1);
        PD.sol_time = [tprofile(1) tprofile(end)];
        PD.init_dist = PD.calc_dist(1);
else    
    
    try
        [PD.calc_time, PD.calc_dist, PD.calc_conc, PD.calc_massmedium] = PBESolver(PD);
    catch ME
        keyboard
        error('ProfileManager:tryconsttemp:PBESolverfail',...
            'PBESolver failed to integrate your problem.')
    end
end