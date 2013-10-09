function [TIME,Y] = hires(PD)
%% [TIME,Y] = hires(PD) High Resolution method for Nucleation and Growth
% Solves the PBE according to a High Resolution method (cf. e.g. Gunawan, R.; Fusman, I.; Braatz, R. D. AIChE Journal 2004, 50, 2738–2749).
% This method does not use a standard ODE solver but rather uses the CFL
% condition + some heuristics to determine the time step size. 
% In general, the high resolution method yields quite accurate results
% while having by far the lowest computational cost. The method is particularly well
% suited for problems with discontinuous distributions (when compared to
% central differences) and virtually never exhibits oscillations.
% Note that, by design, this method has negligible mass balance error,
% however this does not necessarily mean that all results are 100%
% accurate.
% If you are unhappy with the result consider the following options:
% - Increase the number of grid points
% - Decrease reltol {1e-2} and abstol {1e-2} [sol_options] (particularly when encountering oscillations in S!)
% - Set mlim to zero (boundary box finder threshold) {eps} [sol_options] (Important for discontinuous distributions and when mu0 loss is intolerable!)
% - Change the flux limiter function via fluxlim {vanLeer} [sol_options] (see Phifinder for details)
% - Use another method

%% Initial values | Local Variables
% Local time (initial time)
t = PD.sol_time(1);

% Initial time
TIME = t; % is a vector 

% Initial concentration
c = PD.init_conc;  

% Initial temperature
T = PD.Tprofile(TIME);

% Initial solvent + antisolvent mass
m = PD.init_massmedium+PD.ASprofile(TIME)-PD.ASprofile(0);

% Initial mass fraction antisolvent
xm = PD.ASprofile(TIME)/m;

% Initial supersaturation
cs = PD.solubility(T,xm);

% Initial supersaturation
S = c/cs;

% Distribution and Grid
y = PD.init_dist.y(:);
Dy = diff(PD.init_dist.boundaries);
F = PD.init_dist.F; 
Y(1,:) = [F(:)' c]; % Output matrix



% local density function (is padded with zeros)
F_dummy = [0;0;F(:);0]; 
 

%% Tolerances and options
% Default tolerances
mlim = eps; % anteil von partikeln, welcher fuer die simulation vernachlaessigt werden darf --> HAS TO BE 0 WHEN WORKING WITH HEAVYSIDE FUNCTION IN F0
Stol = 1e-2; % tolerance in S (abstol)
ctol = 1e-2; % tolerance in relative change of c and cs (reltol)
fluxlim = 'vanleer';

% if user has set tolerances, use them
if ~isempty(PD.sol_options)
    if ~isempty(find(strcmpi(PD.sol_options,'abstol'),1))
        Stol = PD.sol_options{find(strcmpi(PD.sol_options,'abstol'),1)+1};
    end
    if ~isempty(find(strcmpi(PD.sol_options,'reltol'),1))
        ctol = PD.sol_options{find(strcmpi(PD.sol_options,'reltol'),1)+1};
    end
    if ~isempty(find(strcmpi(PD.sol_options,'mlim'),1))
        mlim = PD.sol_options{find(strcmpi(PD.sol_options,'mlim'),1)+1};
    end
    if ~isempty(find(strcmpi(PD.sol_options,'fluxlim'),1))
        fluxlim = PD.sol_options{find(strcmpi(PD.sol_options,'fluxlim'),1)+1};
    end
end

%% Integration
flagdt = 0; % flag if time step was just rejected (skip time step evaluation)
Dtlast = inf; % last time step

while t<PD.sol_time(end)
       % Growth rate
    if S>=1
        G = PD.growthrate(S,T,PD.init_dist.boundaries(2:end)); % in the high resolution method, the growth rate is evaluated AT THE BOUNDARIES of the bins
    else % dissolution (evaluate at lower boundaries)
        G = PD.growthrate(S,T,PD.init_dist.boundaries(1:end-1)); 
    end

    % Autotimestepsizer based on CFL condition (eq. 24 in Gunawan 2004)
    if flagdt==0
        GI          =   boundingBoxFinder(F_dummy(3:end-1),mlim); % set of indices encapsulating 1-mlim of the distribution
        [~,I]       =   max(abs(G(GI(1):GI(2)-1))./Dy(GI(1):GI(2)-1)); % CFL condition
        Dt          =   abs(Dy(I)/G(I));

        nexttline   =   PD.sol_time(find(PD.sol_time>t,1,'first')); %make sure you hit time points in sol_time vector
        Dttline     =   nexttline-t;

        Dt          =   min([Dtlast*1.5 ... % not more than 3/2 times last time step (acceleration limit)
            max([Dt (PD.sol_time(end)-PD.sol_time(1))*1e-5]) ... % respect CFL, but maximum 1e5 time steps
            Dttline ... % make sure you hit that point
            1/10*(PD.sol_time(end)-PD.sol_time(1))]);  % minimum 10 time steps  
        Dtlast      =   Dt; % save last time step
        F_dummy0      =   F_dummy; % save current distribution
    end %flagdt
    
    % Current mass flow rate antisolvent (evaluated using simplistic FD)
    Q = (PD.ASprofile(t+Dt)-PD.ASprofile(t))/Dt;


    t           =   t+Dt; % update time step

    % mini failsafe
    if PD.sol_time(end)-t     <   1e-12
        t   =   PD.sol_time(end);
    end

    %% Growth     
    if abs(c-cs)>eps 
        F_dummy = hiResGrowth(F_dummy0,G,Dt,Dy,GI,fluxlim); % let it grow, let it grow, let it grow
    else
        F_dummy = F_dummy0;
    end
    
    %% Nucleation
    if  c>cs % nucleation can never occur for S<=1
        if nargin(PD.nucleationrate) == 3 % case where nucleation depends on a moment
            dist = Distribution(y,F_dummy(3:end-1));
            J = PD.nucleationrate(S,T,dist);
        else % nucleation depends only on S and T
            J = PD.nucleationrate(S,T);
        end
        F_dummy(3)  = F_dummy(3) + J/Dy(1)*Dt; 
    end

    %% Calculation of concentration loss due to nucleation and growth
    Deltac = sum((F_dummy(3:end-1)-F_dummy0(3:end-1)).*y.^3.*...
        Dy(:)*PD.rhoc.*PD.kv);
    
    c_dummy    =    c-Deltac;
    
    %% Dilution
    F_dummy = F_dummy *m/(m + Q* Dt); % dilution due to addition of AS
    c_dummy = c_dummy*m/(m+Q*Dt);       
   
    T_dummy    =    PD.Tprofile(t); % next temperature
    xm_dummy   =    PD.ASprofile(t)/m; % next mass fraction antisolvent
    cs_dummy   =    PD.solubility(T_dummy,xm_dummy);     %Solubility
    
    %% Calculate relative and absolute changes    
    DeltaS      =   c_dummy./cs_dummy-c./cs;
    crel        =   abs((c_dummy-c)/c);
    csrel       =   abs((cs_dummy-cs)/cs);
    Qrel        =   Q*Dt/m;

    % Check if result is (superficially) reasonable
    if  (c_dummy>0 && sum(-F_dummy(F_dummy<0))<sum(F_dummy(F_dummy>0))*1e-2 &&...
            abs(DeltaS)<Stol && ((DeltaS>=0 && (crel<ctol && csrel <ctol)) || (DeltaS<0) && (crel<ctol*10 && csrel <ctol)) ...
        && Qrel < Stol...
        ||    flagdt >= 20)

        if flagdt >= 20
            F_dummy(F_dummy<0)=0; % hackebeil
        end

        % Update local variables
        c       =   c_dummy;
        cs      =   cs_dummy;
        T       =   T_dummy;
        S       =   c_dummy/cs_dummy;
        m       =   PD.init_massmedium+PD.ASprofile(t)-PD.ASprofile(0);

        flagdt  =   0; % reset flag

        % save results
        TIME    =   [TIME t]; %#ok<*AGROW>
        F = arrayCrop(F_dummy,[3;length(F_dummy)-1])';
        Y(end+1,:) = [F(:)' c];

    else % results violate tolerances and conditions
        % Use a smaller timestep and repeat everything
        F_dummy =   F_dummy0;
        t       =   t-Dt;
        flagdt  =   flagdt+1;
        Dt      =   Dt/3;
    end


end

end



