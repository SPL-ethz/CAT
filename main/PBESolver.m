function [SolutionTimes SolutionDists SolutionConc] = PBESolver (PD)

% PBESOLVER
%
% Solve PBEs, like a boss

% Setup Solver - currently central difference and moving pivot in separate
% caes of a switch

solvefun = str2func(PD.sol_method);
solvefun = @(t,X) solvefun(t,X,PD);

switch PD.sol_method

    case 'movingpivot'
        dL = 10; % critical bin size for event listener
        X0 = [PD.init_dist.F.*diff(PD.init_dist.boundaries) ...
            PD.init_dist.y PD.init_dist.boundaries PD.init_conc];
        tstart = PD.sol_time(1); % local start time
        tend = PD.sol_time(end); % overall end time

        % if nucleation is present, bins are addded when the first bin
        % becomes too big
        ODEoptions = odeset('RelTol', 1e-8, 'Events',@(t,x) EventAddBin(t,x,dL));
%         keyboard
        nt = length(PD.sol_time);
        if(nt > 2)
            SolutionTimes = zeros(nt,1);
            SolutionConc = zeros(nt,1);
            SolutionDists = repmat(Distribution(),1,nt);
        end
        
        SolutionTimes(1) = tstart;
        SolutionConc(1) = PD.init_conc;
        SolutionDists(1) = PD.init_dist;
        
        X = X0;
        s = 1;
        while tstart<tend 
            ts = PD.sol_time(PD.sol_time > tstart & PD.sol_time < tend);
            if tstart ~= PD.sol_time(1)
                X0 = addBin(X(end,:)'); 
            end
            % Solve until the next event where the nucleation bin becomes to big (as defined by dL)
            [T,X] = ode15s(solvefun, [tstart ts tend],X0, ODEoptions);
            
            if(nt > 2)
                [~, I] = intersect(T,PD.sol_time);
            else
                I = 2:length(T);
            end
            
            tstart = T(end); 
            
            % Break up the output to make it easier to assign. 
            nBins = (size(X,2)-2)/3;            
            y = X(:,nBins+1:2*nBins);            
            boundaries = X(:,2*nBins+1:3*nBins+1);
%             keyboard
            F = X(:,1:nBins)./diff(boundaries,1,2);
            F(isnan(F)) = 0;
            
            for i = 1:length(I)
                SolutionTimes(s+i) = T(I(i));
                SolutionConc(s+i) = X(I(i),3*nBins+2);   
                SolutionDists(s+i) = Distribution( y(I(i),:) , F(I(i),:), boundaries(I(i),:) );
            end % for
            s = s + length(I);
        end %while        
    case 'centraldifference'
        y = PD.init_dist.y;
        X0 = [PD.init_dist.F  PD.init_conc];
        
        [SolutionTimes,X_out] = ode15s(solvefun , PD.sol_time , X0 );
        
        % Create solution        
        SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  %# Pre-Allocation for speed       
        for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution( y , X_out(i,1:length(y)), PD.init_dist.boundaries );
        end % for
        SolutionConc = X_out(:,end);

    case 'hires'

        [SolutionTimes,X_out] = hires(PD);
        
        % Create solution        
        SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  %# Pre-Allocation for speed       
        for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution( PD.init_dist.y, X_out(i,1:length(PD.init_dist.y)),PD.init_dist.boundaries );
        end % for
        SolutionConc = X_out(:,end);

end %switch


end % function