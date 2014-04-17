function [O] = solve(O)
%% CAT.solve
% solve is a method of the CAT class. It integrates the problem according to
% the settings of the CAT object and fills in the calc-properties of the
% same object. The method does not accept any additional inputs, all options 
% should therefore be set via the property fields of the object.
% Solve is structured into two nested parts: An upper layer that handles
% discontinuities in the supplied AS or T profiles and a lower layer (PBEsolver)
% that handles the individual solvers.
%% Solve
O.calc_dist = Distribution;
if ~isempty(O.init_seed)
    O.init_dist.mass = [O.init_seed O.kv O.rhoc O.init_massmedium];
end
sol_time = O.sol_time;
if ~isempty(O.tNodes)
    for i = 2:length(O.tNodes) % make sure you hit the different nodes of the non-smooth profiles
        
        O.sol_time = [O.tNodes(i-1) sol_time(sol_time>O.tNodes(i-1) & sol_time<O.tNodes(i)) O.tNodes(i)];  
        
        try
            [a b c] = PBESolver(O);
        catch ME
            keyboard
            error('solve:tryconsttemp:PBESolverfail',...
            'PBESolver failed to integrate your problem.')
        end
        
        O.calc_time(end+1:end+length(a)) = a;
        O.calc_dist(end+1:end+length(b)) = b;    
        O.calc_conc(end+1:end+length(a)) = c; 

        O.init_dist = O.calc_dist(end);
        O.init_conc = O.calc_conc(end);

    end % for

        O.calc_dist = O.calc_dist(2:end);
        O.sol_time = sol_time;
        O.init_dist = O.calc_dist(1);
else    
    
    try
        PBESolver(O);
    catch ME
        error('solve:tryconsttemp:PBESolverfail',...
            'PBESolver failed to integrate your problem. Message: %s',ME.message)
    end
end

O.init_conc = O.calc_conc(1);
O.init_dist = O.calc_dist(1);

if length(O.sol_time)>2
    [~,I] = intersect(O.calc_time,O.sol_time);
    O.calc_time = O.calc_time(I);
    O.calc_dist = O.calc_dist(I);
    O.calc_conc = O.calc_conc(I);
end

%% Check Mass balance
if any(O.massbal > 5)
   warning('ProfileManager:massbalcheck:largeerror',...
                    'Your mass balance error is unusually large (%4.2f%%). Check validity of equations and consider increasing the number of bins.',max(O.massbal)); 
end

end

function PBESolver (O)

% PBESOLVER
%
% Handles individual solvers and returns (local) results

solvefun = str2func(['@(t,X) ' O.sol_method '(t,X,O)']);

switch O.sol_method

    case 'movingpivot'
        O.movingpivot;       
        
    case 'centraldifference'
        
        options = O.sol_options;
        if isempty(O.sol_options)
            options = odeset(options,'reltol',1e-6);
        end
        
        X0 = [O.init_dist.F, O.init_conc];
        
        [SolutionTimes,X_out] = ode15s(solvefun , O.sol_time , X0 ,options);

    case 'hires'
        
        [SolutionTimes,X_out] = hires(O);
         
end %switch

% For CD and HR transform result-arrays into appropriate output structure
if ~strcmpi(O.sol_method,'movingpivot')
    SolutionConc = X_out(:,end);
    SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  % Pre-Allocation for speed               
    for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution( O.init_dist.y, X_out(i,1:length(O.init_dist.y)),O.init_dist.boundaries );
    end % for
end

end % function