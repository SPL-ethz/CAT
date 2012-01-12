function [SolutionTimes SolutionDists SolutionConc SolutionTemp SolutionVolume] = PBESolver (PD)

% PBESOLVER
%
% Solve PBEs, like a boss

% Setup Solver - currently central difference and moving pivot in separate
% caes of a switch

solvefun = str2func(PD.sol_method);
solvefun = @(t,X) solvefun(t,X,PD);

switch PD.sol_method
    case 'hires'
        finput.exp.f0 = PD.init_dist.F;

        if ~isempty(PD.init_dist.boundaries)
            finput.num.boundaries.dim1  = PD.init_dist.boundaries;
            finput.num.ngrid            = length(PD.init_dist.boundaries)-1;
        elseif ~isempty(PD.init_dist.y)
            finput.num.y.dim1  = PD.init_dist.y;
            finput.num.ngrid   = length(PD.init_dist.y);
        end

        output=HRFL_1dim(finput,PD);
        output.PSD
        SolutionTimes = output.time;
        SolutionConc = output.c;
        
        for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution(output.PSD.xp.dim1,output.PSD.F(:,i),output.PSD.xb.dim1);
        end
 
    case 'movingpivot'
        dL = 20e-6;
        X0 = [PD.init_dist.F PD.init_dist.y PD.init_dist.boundaries  PD.init_conc PD.init_temp PD.init_volume];
        tstart = PD.sol_time(1);
        tend = PD.sol_time(end);
        
        ODEoptions = odeset;
        
        % if nucleation is present, bins are addded when the first bin
        % becomes too big
        if PD.nucleationrate(PD.init_conc,PD.init_temp) > 0
            ODEoptions = odeset('Events',@(t,x) addBinEvent(t,x,dL));
            X0 = addBin(X0(:)); %since nucleation is present, add an empty bin
        end
        
        nt = length(PD.sol_time);
        if(nt > 2)
            SolutionTimes = zeros(nt,1);
            SolutionConc = zeros(nt,1);
            SolutionTemp = zeros(nt,1);
            SolutionVolume = zeros(nt,1);
            SolutionDists = repmat(Distribution(),1,nt);
        end
        
        SolutionTimes(1) = tstart;
        SolutionConc(1) = PD.init_conc;
        SolutionDists(1) = PD.init_dist;
        SolutionTemp(1) = PD.init_temp;
        SolutionVolume(1) = PD.init_volume;
        
        s = 1;
        while tstart<tend 
            ts = PD.sol_time(PD.sol_time > tstart & PD.sol_time < tend);
            
            % Solve until the next event where the nucleation bin becomes to big (as defined by parameters.dL)
            [T,X] = ode15s(solvefun, [tstart ts tend],X0, ODEoptions);
            
            if(nt > 2)
                [~, I] = intersect(T,PD.sol_time);
            else
                I = 1:length(T);
            end
            
            tstart = T(end); 
            
            % Break up the output to make it easier to assign. 
            nBins = (size(X,2)-4)/3;
            F = X(:,1:nBins);
            y = X(:,nBins+1:2*nBins);            
            boundaries = X(:,2*nBins+1:3*nBins+1);    
            
            for i = 1:length(I)
                SolutionTimes(s+i) = T(I(i));
                SolutionConc(s+i) = X(I(i),3*nBins+2);   
                SolutionTemp(s+i) = X(I(i),3*nBins+3);   
                SolutionVolume(s+i) = X(I(i),3*nBins+4);   
                SolutionDists(s+i) = Distribution( y(I(i),:) , F(I(i),:), boundaries(I(i),:) );
            end % for
            
            ODEoptions = odeset(ODEoptions,'InitialStep',T(end)-T(end-1)); %choosing step size based on the last step size
            X0 = addBin(X(end,:)'); 
            s = s + length(I);
        end %while        
    case 'centraldifference'
        y = PD.init_dist.y;
        X0 = [PD.init_dist.F  PD.init_conc];
        
        [SolutionTimes,X_out] = ode15s(solvefun , PD.sol_time , X0 );
        
        % Create solution        
        SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  %# Pre-Allocation for speed       
        for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution( y , X_out(i,1:length(y)) );
        end % for
        SolutionConc = X_out(:,end);
        SolutionTemp = ones(size(SolutionConc))*PD.init_temp;
        SolutionVolume = ones(size(SolutionConc))*PD.init_volume;
end %switch
end % function