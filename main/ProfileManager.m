function [PD] = ProfileManager(PD)

flagTp = 0;flagQp = 0;tprofile=[];
%% Solve
PD.calc_dist = Distribution;
if ~isempty(PD.Tprofile) && ~isempty(PD.tTprofile) && length(PD.Tprofile) == length(PD.tTprofile)
    flagTp = 1;
elseif ~isempty(PD.coolingrateprofile) && ~isempty(PD.tTprofile) && length(PD.coolingrateprofile) == length(PD.tTprofile)-1
    flagTp = 2;
end

if ~isempty(PD.ASadditionrateprofile) && ~isempty(PD.tQprofile) && length(PD.ASadditionrateprofile) == length(PD.tQprofile)-1
    flagQp = 1;
end

if flagTp ~=0 && flagQp ~= 0
    if PD.tTprofile(1) ~= PD.tQprofile(1) || PD.tTprofile(end) ~= PD.tQprofile(end)
        warning('ProfileManager:tTQprofile:inconsistenttprofiles',...
            'The Temperature and AS Addition Rate profiles are inconsistent');
    end
end
        
if flagTp ~= 0 || flagQp ~= 0
    tprofile = sort(unique([PD.tTprofile PD.tQprofile]),'ascend');

    for i = 2:length(tprofile);
        
        if flagTp == 1
            iT = find(PD.tTprofile<tprofile(i),1,'last');
            PD.coolingrate = (PD.Tprofile(iT+1)-PD.Tprofile(iT))/(PD.tTprofile(iT+1)-PD.tTprofile(iT));
        elseif flagTp == 2
            iT = find(PD.tTprofile<tprofile(i),1,'last');
            PD.coolingrate = PD.coolingrateprofile{iT};
        end
        
        if flagQp == 1
            iQ = find(PD.tQprofile<tprofile(i),1,'last');
            PD.ASadditionrate = PD.ASadditionrateprofile{iQ};
        end
        
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
    end % for

        PD.calc_dist = PD.calc_dist(2:end);
        PD.init_temp = PD.calc_temp(1);
        PD.sol_time = [tprofile(1) tprofile(end)];
        PD.init_dist = PD.calc_dist(1);
else
    
    if xor(isempty(PD.Tprofile),isempty(PD.tTprofile)) || xor(isempty(PD.coolingrateprofile),isempty(PD.tTprofile))
        warning('ProfileManager:profile:emptyprofile',...
            'Either T, coolingrate or t profile is empty. \n Attempting to use constant Temperature')
        
    elseif ~isempty(PD.Tprofile) && ~isempty(PD.tTprofile) && length(PD.Tprofile) ~= length(PD.tTprofile)
        warning('ProfileManager:Ttprofiles:unequallength',...
            'T and t profiles have unequal size. \n Attempting to use constant Temperature')
        
    elseif ~isempty(PD.coolingrateprofile) && ~isempty(PD.tTprofile) && length(PD.coolingrateprofile) ~= length(PD.tTprofile)-1
    warning('ProfileManager:ctprofiles:unequallength',...
        'Coolingrate and t profiles have inappropriate sizes (cooling rate must be length(tTprofile)-1. \n Attempting to use constant Temperature')
    end
        
    
    try
        [PD.calc_time, PD.calc_dist, PD.calc_conc, PD.calc_temp, PD.calc_volume] = PBESolver(PD);
    catch
        error('ProfileManager:tryconsttemp:PBESolverfail',...
            'PBESolver failed to integrate your problem.')
    end
end