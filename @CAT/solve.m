function [O] = solve(O)

%% Solve
O.calc_dist = Distribution;
O.init_dist.mass = [O.init_seed O.kv O.rhoc O.init_massmedium];
sol_time = O.sol_time;
if ~isempty(O.tNodes)
    for i = 2:length(O.tNodes) % make sure you hit the different nodes of the non-smooth profiles
        
        O.sol_time = [O.tNodes(i-1) sol_time(sol_time>O.tNodes(i-1) & sol_time<O.tNodes(i)) O.tNodes(i)];  
        
        [a b c] = PBESolver(O);
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
        
        [O.calc_time, O.calc_dist, O.calc_conc] = PBESolver(O);
    catch ME
        keyboard
        error('ProfileManager:tryconsttemp:PBESolverfail',...
            'PBESolver failed to integrate your problem.')
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

function [SolutionTimes SolutionDists SolutionConc] = PBESolver (O)

% PBESOLVER
%
% Solve PBEs, like a boss

solvefun = str2func(O.sol_method);
solvefun = @(t,X) solvefun(t,X,O);

switch O.sol_method

    case 'movingpivot'
        if ~isempty(O.sol_options) && ~isempty(find(strcmpi(O.sol_options,'dL'),1))
            dL = O.sol_options(find(strcmpi(O.sol_options,'dL'),1)+1);
        else
            dL = 10; % critical bin size for event listener
        end
        X0 = [O.init_dist.F.*diff(O.init_dist.boundaries) ...
            O.init_dist.y O.init_dist.boundaries O.init_conc];
        tstart = O.sol_time(1); % local start time
        tend = O.sol_time(end); % overall end time

        % if nucleation is present, bins are addded when the first bin
        % becomes too big
        options = O.sol_options;
        if isempty(O.sol_options)
            options = odeset(options,'Events',@(t,x) EventBin(t,x,dL),'reltol',1e-6);            
        else
            options = odeset(options,'Events',@(t,x) EventBin(t,x,dL));
        end

        
        SolutionTimes = []; SolutionConc = [];
        s=0;
        while tstart<tend 
            ts = O.sol_time(O.sol_time > tstart & O.sol_time < tend);
            
            % Solve until the next event where the nucleation bin becomes to big (as defined by dL)
            [TIME,X_out] = ode15s(solvefun, [tstart ts tend],X0, options);
            
            X_out(X_out<0) = 0;
            
            nBins = (size(X_out,2)-2)/3;
            
            if X_out(end,2*nBins+1)>0
                X0 = addBin(X_out(end,:)'); 
            elseif X_out(end,nBins+1)<=0
                X0 = removeBin(X_out(end,:)'); 
            end
                                   
            F = X_out(:,1:nBins)./diff(X_out(:,2*nBins+1:3*nBins+1),1,2); F(isnan(F)) = 0;
            
            SolutionTimes = [SolutionTimes;TIME(:)]; %#ok<AGROW>
            SolutionConc = [SolutionConc;X_out(:,end)]; %#ok<AGROW>
            for i = 1:length(TIME)
                SolutionDists(s+i) = Distribution( X_out(i,nBins+1:2*nBins),...
                    F(i,:),...
                    X_out(i,2*nBins+1:3*nBins+1) ); %#ok<AGROW>
            end
            s = s+length(TIME);
            tstart = TIME(end); 
            
        end %while        
        
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

if ~strcmpi(O.sol_method,'movingpivot')
    SolutionConc = X_out(:,end);
    SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  %# Pre-Allocation for speed               
    for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution( O.init_dist.y, X_out(i,1:length(O.init_dist.y)),O.init_dist.boundaries );
    end % for
end

end % function