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

        output = hires(finput,PD);
        SolutionTimes = output.time;
        SolutionConc = output.c;
        SolutionTemp = output.Temp;
        SolutionVolume = output.Volume;
        
        SolutionDists = repmat(Distribution(),1,length(SolutionTimes));  %# Pre-Allocation for speed 
        for i = 1:length(SolutionTimes)
            SolutionDists(i) = Distribution(output.PSD.xp.dim1,output.PSD.F(:,i),output.PSD.xb.dim1);
        end
 
    case 'movingpivot'
        dL = 50e-6;
        X0 = [PD.init_dist.F.*(PD.init_dist.boundaries(2:end)-PD.init_dist.boundaries(1:end-1)) ...
            PD.init_dist.y PD.init_dist.boundaries  PD.init_conc PD.init_temp PD.init_volume];
        tstart = PD.sol_time(1);
        tend = PD.sol_time(end);
        
        
        
        % if nucleation is present, bins are addded when the first bin
        % becomes too big
        ODEoptions = odeset('RelTol', 1e-8);         
%         if PD.nucleationrate(PD.init_conc,PD.init_temp) > 0
            ODEoptions = odeset(ODEoptions, 'Events',@(t,x) EventAddBin(t,x,dL));
%        end

        
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
        
        X = X0;
        s = 1;
        while tstart<tend 
            ts = PD.sol_time(PD.sol_time > tstart & PD.sol_time < tend);
            X0 = addBin(X(end,:)'); 
            % Solve until the next event where the nucleation bin becomes to big (as defined by dL)
            [T,X] = ode15s(solvefun, [tstart ts tend],X0, ODEoptions);
            
            if(nt > 2)
                [~, I] = intersect(T,PD.sol_time);
            else
                I = 2:length(T);
            end
            
            tstart = T(end); 
            
            % Break up the output to make it easier to assign. 
            nBins = (size(X,2)-4)/3;            
            y = X(:,nBins+1:2*nBins);            
            boundaries = X(:,2*nBins+1:3*nBins+1);
            F = X(:,1:nBins)./(boundaries(:,2:end)-boundaries(:,1:end-1));
            F(isnan(F)) = 0;
            
            for i = 1:length(I)
                SolutionTimes(s+i) = T(I(i));
                SolutionConc(s+i) = X(I(i),3*nBins+2);   
                SolutionTemp(s+i) = X(I(i),3*nBins+3);   
                SolutionVolume(s+i) = X(I(i),3*nBins+4);   
                SolutionDists(s+i) = Distribution( y(I(i),:) , F(I(i),:), boundaries(I(i),:) );
            end % for
%             ODEoptions = odeset(ODEoptions,'InitialStep',T(nt)-T(end-1),...
%                            'MaxStep',T(end)-T(1));
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