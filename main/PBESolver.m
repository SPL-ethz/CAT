function [SolutionTimes SolutionDists SolutionConc] = PBESolver (PD)

% PBESOLVER
%
% Solve PBEs, like a boss

solvefun = str2func(PD.sol_method);
solvefun = @(t,X) solvefun(t,X,PD);

switch PD.sol_method

    case 'movingpivot'
        if ~isempty(PD.sol_options) && ~isempty(find(strcmpi(PD.sol_options,'dL'),1))
            dL = PD.sol_options(find(strcmpi(PD.sol_options,'dL'),1)+1);
        else
            dL = 10; % critical bin size for event listener
        end
        X0 = [PD.init_dist.F.*diff(PD.init_dist.boundaries) ...
            PD.init_dist.y PD.init_dist.boundaries PD.init_conc];
        tstart = PD.sol_time(1); % local start time
        tend = PD.sol_time(end); % overall end time

        % if nucleation is present, bins are addded when the first bin
        % becomes too big
        options = PD.sol_options;
        if isempty(PD.sol_options)
            options = odeset(options,'Events',@(t,x) EventBin(t,x,dL),'reltol',1e-6);            
        else
            options = odeset(options,'Events',@(t,x) EventBin(t,x,dL));
        end

        
        SolutionTimes = []; SolutionConc = [];
        s=0;
        while tstart<tend 
            ts = PD.sol_time(PD.sol_time > tstart & PD.sol_time < tend);
            
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
        
        options = PD.sol_options;
        if isempty(PD.sol_options)
            options = odeset(options,'reltol',1e-6);
        end
        
        X0 = [PD.init_dist.F, PD.init_conc];
        
        [SolutionTimes,X_out] = ode15s(solvefun , PD.sol_time , X0 ,options);

    case 'hires'
        
        [SolutionTimes,X_out] = hires(PD);
         
end %switch

if ~strcmpi(PD.sol_method,'movingpivot')
    SolutionConc = X_out(:,end);
    SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  %# Pre-Allocation for speed               
    for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution( PD.init_dist.y, X_out(i,1:length(PD.init_dist.y)),PD.init_dist.boundaries );
    end % for
end

end % function