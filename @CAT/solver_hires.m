function [mbflag] = solver_hires(O)
%% solver_hires(O)
% High Resolution method for Nucleation and Growth
% Solves the PBE according to a High Resolution method (cf. e.g. Gunawan, R.; Fusman, I.; Braatz, R. D. AIChE Journal 2004, 50, 2738�2749).
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
t = O.sol_time(1);

% Initial time
TIME = t; % is a vector 

% Initial concentration
c = O.init_conc;  

% Initial temperature
T = O.Tprofile(TIME);

% Initial solvent + antisolvent mass
m = O.init_massmedium+O.ASprofile(TIME)-O.ASprofile(0);

% Initial mass fraction antisolvent
xm = O.ASprofile(TIME)/m;

% Initial supersaturation
cs = O.solubility(T,xm);

% Initial supersaturation
S = c/cs;

% Distribution and Grid
y = O.init_dist.y(:);
Dy = diff(O.init_dist.boundaries);
F = O.init_dist.F; 
Y(1,:) = [F(:)' c]; % Output matrix



% local density function (is padded with zeros)
F_dummy = [0;0;F(:);0]; 
 

%% Tolerances and options
% Default tolerances
mlim = 0; % anteil von partikeln, welcher fuer die simulation vernachlaessigt werden darf --> HAS TO BE 0 WHEN WORKING WITH HEAVYSIDE FUNCTION IN F0
Stol = 1e-2; % tolerance in S (abstol)
ctol = 1e-2; % tolerance in relative change of c and cs (reltol)
fluxlim = 'vanleer';

% if user has set tolerances, use them
if ~isempty(O.sol_options)
    if ~isempty(find(strcmpi(O.sol_options,'abstol'),1))
        Stol = O.sol_options{find(strcmpi(O.sol_options,'abstol'),1)+1};
    end
    if ~isempty(find(strcmpi(O.sol_options,'reltol'),1))
        ctol = O.sol_options{find(strcmpi(O.sol_options,'reltol'),1)+1};
    end
    if ~isempty(find(strcmpi(O.sol_options,'mlim'),1))
        mlim = O.sol_options{find(strcmpi(O.sol_options,'mlim'),1)+1};
    end
    if ~isempty(find(strcmpi(O.sol_options,'fluxlim'),1))
        fluxlim = O.sol_options{find(strcmpi(O.sol_options,'fluxlim'),1)+1};
    end
end

%% Integration
flagdt = 0; % flag if time step was just rejected (skip time step evaluation)
Dtlast = inf; % last time step
if isempty(O.tNodes)
    tFinal = O.sol_time(end);
else
    tFinal = O.tNodes(end);
end

while t<O.sol_time(end)
       % Growth rate
    if S>=1
        G = O.growthrate(S,T,O.init_dist.boundaries(2:end)); % in the high resolution method, the growth rate is evaluated AT THE BOUNDARIES of the bins
    else % dissolution (evaluate at lower boundaries)
        G = O.growthrate(S,T,O.init_dist.boundaries(1:end-1)); 
    end

    % Autotimestepsizer based on CFL condition (eq. 24 in Gunawan 2004)
    if flagdt==0
        GI          =   boundingBoxFinder(F_dummy(3:end-1),mlim); % set of indices encapsulating 1-mlim of the distribution
        [~,I]       =   max(abs(G(GI(1):GI(2)-1))./Dy(GI(1):GI(2)-1)); % CFL condition
        Dt          =   abs(Dy(I)/G(I));

        nexttline   =   O.sol_time(find(O.sol_time>t,1,'first')); %make sure you hit time points in sol_time vector
        Dttline     =   nexttline-t;

        Dt          =   min([Dtlast*1.5 ... % not more than 3/2 times last time step (acceleration limit)
            max([Dt (O.sol_time(end)-O.sol_time(1))*1e-5]) ... % respect CFL, but maximum 1e5 time steps
            Dttline ... % make sure you hit that point
            1/10*(O.sol_time(end)-O.sol_time(1))]);  % minimum 10 time steps  
        Dtlast      =   Dt; % save last time step
        F_dummy0      =   F_dummy; % save current distribution
    end %flagdt
    
    % Current mass flow rate antisolvent (evaluated using simplistic FD)
    Q = (O.ASprofile(t+Dt)-O.ASprofile(t))/Dt;


    t           =   t+Dt; % update time step

    % mini failsafe
    if O.sol_time(end)-t     <   1e-12
        t   =   O.sol_time(end);
    end

    %% Growth     
    if abs(c-cs)>eps 
        F_dummy = hiResGrowth(F_dummy0,G,Dt,Dy,GI,fluxlim); % let it grow, let it grow, let it grow
    else
        F_dummy = F_dummy0;
    end
    
    %% Nucleation
    if  c>cs % nucleation can never occur for S<=1
        if nargin(O.nucleationrate) > 2 % case where nucleation depends on a moment
            dist = Distribution(y,F_dummy(3:end-1));
            J = O.nucleationrate(S,T,dist);
        else % nucleation depends only on S and T
            J = O.nucleationrate(S,T);
        end
        F_dummy(3)  = F_dummy(3) + J/Dy(1)*Dt; 
    end

    %% Calculation of concentration loss due to nucleation and growth
    Deltac = sum((F_dummy(3:end-1)-F_dummy0(3:end-1)).*y.^3.*...
        Dy(:)*O.rhoc.*O.kv);
    
    c_dummy    =    c-Deltac;
    
    %% Dilution
    F_dummy = F_dummy *m/(m + Q* Dt); % dilution due to addition of AS
    c_dummy = c_dummy*m/(m+Q*Dt);       
   
    T_dummy    =    O.Tprofile(t); % next temperature
    xm_dummy   =    O.ASprofile(t)/m; % next mass fraction antisolvent
    cs_dummy   =    O.solubility(T_dummy,xm_dummy);     %Solubility
    
    %% Calculate relative and absolute changes    
    DeltaS      =   c_dummy./cs_dummy-c./cs;
    crel        =   abs((c_dummy-c)/c);
    csrel       =   abs((cs_dummy-cs)/cs);
    Qrel        =   Q*Dt/m;

    % Check if result is (superficially) reasonable
    if  (c_dummy>0 && sum(-F_dummy(F_dummy<0))<=sum(F_dummy(F_dummy>0))*1e-2 &&...
            sum(F_dummy)==0 || (abs(DeltaS)<Stol && ((DeltaS>=0 && (crel<ctol && csrel <ctol)) || (DeltaS<0) && (crel<ctol*10 && csrel <ctol)) ...
        && Qrel < Stol)...
        ||    flagdt >= 20)

        if flagdt >= 20
            F_dummy(F_dummy<0)=0; % hackebeil
        end

        % Update local variables
        c       =   c_dummy;
        cs      =   cs_dummy;
        T       =   T_dummy;
        S       =   c_dummy/cs_dummy;
        m       =   O.init_massmedium+O.ASprofile(t)-O.ASprofile(0);

        flagdt  =   0; % reset flag

        % save results
        TIME    =   [TIME t]; %#ok<*AGROW>
        F = arrayCrop(F_dummy,[3;length(F_dummy)-1])';
        Y(end+1,:) = [F(:)' c];
        
        if findall(0,'name','Looking at CATs')

            if floor(t/tFinal/0.05)>str2num(get(gca,'tag'))
                fill([0 t/tFinal t/tFinal 0],[0 0 1 1],'c','edgecolor','none')
                delete(findall(gcf,'type','text'))
                text(0.44,0.5,[num2str(floor(t/tFinal*100),'%2d'),'%'])
                set(gca,'tag',num2str(floor(t/tFinal/0.05)))
                drawnow
            end
        end

    else % results violate tolerances and conditions
        % Use a smaller timestep and repeat everything
        F_dummy =   F_dummy0;
        t       =   t-Dt;
        flagdt  =   flagdt+1;
        Dt      =   Dt/3;
    end

mbflag = 0;

end

%% Assign output

% Transform result-arrays into appropriate output structure
SolutionConc = Y(:,end);
SolutionDists = repmat(Distribution(),1,length(TIME));  % Pre-Allocation for speed
for i = 1:length(TIME)
    SolutionDists(i) = Distribution( O.init_dist.y, Y(i,1:length(O.init_dist.y)),O.init_dist.boundaries );
end % for

O.calc_time = TIME;
O.calc_dist = SolutionDists;
O.calc_conc = SolutionConc;



end

%% Function hiResGrowth

function[Ffull] = hiResGrowth(Ffull,G,Dt,Dy,GI,fluxlim)
   
FI = zeros(size(GI));FI(1) = max([GI(1)-3 1]);FI(2) = min([GI(2)+3 length(Ffull)]); % boundary box for (padded) distribution
Dy = [eps;eps;Dy(:); eps];Dy = Dy(FI(1):FI(2));

F = Ffull(FI(1):FI(2)); F = F(:);

% Find Theta
if any(G>=  0) % growth
    Theta   =   arrayDivision(diff(F,1,1),1,1);
else % dissolution
    F(1:2) = repmat(F(3),[2 1]); % outflow boundary condition p 131 ff leveque
    Theta   =   arrayDivision(diff(F,1,1),1,2);
    Theta   =   circshift(Theta,[-1 0]);
end

Theta(isinf(Theta))     =   0;
Theta(isnan(Theta))     =   2;

Phi                     =   Phifinder(Theta,fluxlim);    % Flux limiter
Phi(isnan(Phi))         =   0;

% Propagate Distribution 
if length(unique(G))==1 % size independent growth
    Gloc = G(1);
    if min(G)>=0

        F(3:end-1) = F(3:end-1)-Dt.*Gloc./Dy(3:end-1).*(F(3:end-1)-F(2:end-2))-Dt.*Gloc./(2*Dy(3:end-1)).*(1-Dt.*Gloc./Dy(3:end-1)).*...
            ((F(4:end)-F(3:end-1)).*Phi(2:end)-(F(3:end-1)-F(2:end-2,:,:,:)).*Phi(1:end-1));

    else
        F(3:end-1) = F(3:end-1)-Dt*Gloc./Dy(4:end).*(F(4:end)-F(3:end-1))-Dt*Gloc./(2*Dy(4:end)).*(1+Dt*Gloc./Dy(4:end)).*...
            ((F(3:end-1)-F(2:end-2)).*Phi(1:end-1)-(F(4:end)-F(3:end-1)).*Phi(2:end));     
    end

elseif length(unique(G))>1 % size depenent growth

    Gloc = [0;0;G(:); 0];Gloc = Gloc(FI(1):FI(2));Gloc=Gloc(:);

    if min(Gloc(:))>=0
        F(3:end-1)=F(3:end-1)-Dt./Dy(3:end-1).*(Gloc(3:end-1).*F(3:end-1)-Gloc(2:end-2).*F(2:end-2))-...
                        (Dt./(2*Dy(3:end-1)).*Gloc(3:end-1).*(1-(Dt./Dy(3:end-1).*Gloc(3:end-1))).*(F(4:end)-F(3:end-1)).*Phi(2:end)-...
                        Dt./(2*Dy(2:end-2)).*Gloc(2:end-2).*(1-Dt./Dy(2:end-2).*Gloc(2:end-2)).*(F(3:end-1)-F(2:end-2)).*Phi(1:end-1));

    else
%         try
        F(3:end-1)=F(3:end-1)-Dt./Dy(4:end).*(Gloc(4:end).*F(4:end)-Gloc(3:end-1).*F(3:end-1))+...          
                    (Dt./(2*Dy(3:end-1)).*Gloc(3:end-1).*(1+(Dt./Dy(3:end-1).*Gloc(3:end-1))).*(F(4:end)-F(3:end-1)).*Phi(2:end)-...
                    Dt./(2*Dy(2:end-2)).*Gloc(2:end-2).*(1+Dt./Dy(2:end-2).*Gloc(2:end-2)).*(F(3:end-1)-F(2:end-2)).*Phi(1:end-1));
%         catch
%             keyboard
%         end
    end

end


Ffull(FI(1):FI(2)) = F;

end % function hiResGrowth

%% Function arrayCrop

function y = arrayCrop(x,cropSize)

sfx = size(x);sfx(sfx==1) = [];
nDim = length(sfx);

ind = {};

for k = 1:nDim
   ind = [ind cropSize(1,k):cropSize(2,k)];
end

y = x(ind{:});

end % function

%% Function arrayDivision

function [Astar] = arrayDivision(A,dim,order)
% divides array A element by element along dimension dim, so that 
% arrayDivision(A,1) gives A(2,:,:)./A(1,:,:) etc. for a 3D array    

%shift the dimensions of the matrix in such a way that the dimension to
%be divised is the last one
if length(size(A))>=3

    A = shiftdim(A,dim);

    %calculate the indeces of elements to be divided in the different
    %layers. variable blub contains the pairs of indices to be divided in
    %two adjacent columns
    s = size(A);
    inc = prod(s(1:end-1));
    basis = linspace(1,inc,inc)';
    increments = (1:inc:numel(A))-1;
    blub = basis*ones(1,length(increments))+ones(inc,1)*increments;

    %do the division
    if order==2
        Astar = A(blub(:,2:end))./A(blub(:,1:end-1));
    elseif order==1
        Astar = A(blub(:,1:end-1))./A(blub(:,2:end));
    end

    %this is the new size of the divided matrix
    s(end) = s(end)-1;

    %reshaping
    Astar = reshape(Astar, s);

    %re-shifting to reflect the initial dimensions

    Astar = shiftdim(Astar,length(s)-dim);

elseif length(size(A))<3
    % do the same also for matrices and vectors
    A=shiftdim(A,dim-1);
    if order == 2
        Astar=A(2:end,:)./A(1:end-1,:);
    elseif order == 1
        Astar=A(1:end-1,:)./A(2:end,:);
    end

    Astar=shiftdim(Astar,dim-1);
end

end % function

%% Function boundingBoxFinder

function [GI] = boundingBoxFinder(F,mlim)

Fcum = cumsum(abs(F));

a0 = find(Fcum./Fcum(end)>mlim,1,'first');
a0(isempty(a0)) = 1; 

af = find(1-Fcum./Fcum(end)<mlim,1,'first');
af(isempty(af)) = length(F(3:end-1,1,1));

GI(1,1) = min([a0 max([1 af-1])]);
GI(1,2) = max([2 af]);

end % function

%% Function Phifinder

function [Phi] = Phifinder(Theta,fluxlim)

% Phifinder(Theta,fluxlim)
% Given a theta, find the corresponding flux limiter. The user can choose
% between 5 different functions for Phi (cf. Leveque 2002). Default is
% VanLeer's function for the limiter
% Theta = Ratio of differences between 3 cells
% fluxlim = String determining function that should be used

switch fluxlim
    case {'vanleer'}
        Phi     =   (abs(Theta)+Theta)./(1+abs(Theta));
    case{'koren'}
        % Check expression here again, might be a little off (different
        % description of theta?)
        Phi = zeros(length(Theta),1);
        for i=1:length(Theta)
            Phi(i)  =   max([0 min([2.*Theta(i) min([1/3+2/3*Theta(i) 2])])]);
        end
        
    case{'minmod'}
        Phi = zeros(length(Theta),1);
        for i=1:length(Theta)
            Phi(i)  =   max([1 Theta(i)]);
        end
        
    case{'superbee'}
        Phi = zeros(length(Theta),1);
        for i=1:length(Theta)
            Phi(i)  =   max([0 min([1 2*Theta(i)]) min([2 Theta(i)])]);
        end
        
    case{'MC'}
        Phi = zeros(length(Theta),1);
        for i=1:length(Theta)
            Phi(i)  =   max([0 min([2 2.*Theta(i) (1+Theta(i))/2])]);
        end
        
    otherwise
        % Default is VanLeer
        Phi=(abs(Theta)+Theta)./(1+abs(Theta));
        
end

end

%%