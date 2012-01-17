function [PD] = ProfileManager(PD)

%% Solve
PD.calc_dist = Distribution;
if ~isempty(PD.Tprofile) && ~isempty(PD.tprofile)
    
    for i = 2:length(PD.Tprofile)
        i
        PD.coolingrate = (PD.Tprofile(i)-PD.Tprofile(i-1))/(PD.tprofile(i)-PD.tprofile(i-1));
        PD.sol_time = [PD.tprofile(i-1) PD.tprofile(i)];  
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
    PD.init_temp = PD.Tprofile(1);
    PD.sol_time = [PD.tprofile(1) PD.tprofile(end)];
    PD.init_dist = PD.calc_dist(1);
    
else
    if xor(isempty(PD.Tprofile),isempty(PD.tprofile))
        warning('ProfileManager:Tprofile:emptyprofile',...
            'Either T or t profile is empty. \n Attempting to use constant Temperature')
    end
        
        try
            [PD.calc_time, PD.calc_dist, PD.calc_conc, PD.calc_temp, PD.calc_volume] = PBESolver(PD);
        catch
            error('ProfileManager:tryconsttemp:PBESolverfail',...
                'PBESolver failed to integrate your problem.')
        end
end