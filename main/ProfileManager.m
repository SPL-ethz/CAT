function [PD] = ProfileManager(PD)

%% Solve
PD.calc_dist = Distribution;
% keyboard
PD.init_dist.mass = [PD.init_seed PD.kv PD.rhoc PD.init_massmedium];

if ~isempty(PD.tNodes)
    for i = 2:length(PD.tNodes) % make sure you hit the different nodes of the non-smooth profiles
        
        PD.sol_time = [PD.tNodes(i-1) PD.tNodes(i)];  
        
        [a b c] = PBESolver(PD);
%         keyboard
        PD.calc_time(end+1:end+length(a)) = a;
        PD.calc_dist(end+1:end+length(b)) = b;    
        PD.calc_conc(end+1:end+length(a)) = c; 

        PD.init_dist = PD.calc_dist(end);
        PD.init_conc = PD.calc_conc(end);

    end % for

        PD.calc_dist = PD.calc_dist(2:end);
        PD.sol_time = [PD.tNodes(1) PD.tNodes(end)];
        PD.init_dist = PD.calc_dist(1);
else    
    
    try
        
        [PD.calc_time, PD.calc_dist, PD.calc_conc] = PBESolver(PD);
    catch ME
        keyboard
        error('ProfileManager:tryconsttemp:PBESolverfail',...
            'PBESolver failed to integrate your problem.')
    end
end

%% Make sanity check of results
% keyboard
if any(PD.massbal > 5)
   warning('ProfileManager:massbalcheck:largeerror',...
                    'Your mass balance error is unusually large (%4.2f%%). Check validity of equations and consider increasing the number of bins.',max(PD.massbal)); 
end