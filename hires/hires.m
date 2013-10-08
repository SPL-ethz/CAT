function [TIME,Y] = hiResSlave(PD,varargin)

%% Setup and Preparation

c = PD.init_conc;  
T = PD.Tprofile(PD.sol_time(1));


TIME = PD.sol_time(1);
%% Integration
flagdt=0;
t=PD.sol_time(1);tcount=1;

fstar = PD.init_dist.F(:);
y = PD.init_dist.y(:);
fstar = [0;0;fstar;0];
sfx = length(fstar);
F = arrayCrop(fstar,[3;sfx-1]);
Y(1,:) = [F(:)' c];
Dy = diff(PD.init_dist.boundaries);


mlim = 0; % anteil von partikeln, welcher fuer die simulation vernachlaessigt werden darf --> HAS TO BE 0 WHEN WORKING WITH HEAVYSIDE FUNCTION IN F0
Stol = 1e-2;
ctol = 1e-2;
if ~isempty(PD.ODEoptions)
    if ~isempty(find(strcmpi(PD.ODEoptions),'abstol'))
        I=find(strcmpi(PD.ODEoptions),'abstol');
        Stol = PD.ODEoptions{I+1};
    end
    if ~isempty(find(strcmpi(PD.ODEoptions),'reltol'))
        I=find(strcmpi(PD.ODEoptions),'reltol');
        ctol = PD.ODEoptions{I+1};
    end
    if ~isempty(find(strcmpi(PD.ODEoptions),'mlim'))
        I=find(strcmpi(PD.ODEoptions),'mlim');
        mlim = PD.ODEoptions{I+1};
    end
end

Dtlast = inf;
while t<PD.sol_time(end)
    % Current solvent + antisolvent mass
    m = PD.init_massmedium+(PD.ASprofile(t)-PD.ASprofile(0));

    % Current mass fraction antisolvent
    xm = PD.ASprofile(t)/m;


    cs = PD.solubility(T,xm);

    % Current supersaturation
    S = c/cs;

    % Current mass flow rate antisolvent (evaluated using simplistic FD)
    Q = (PD.ASprofile(t+1e-6)-PD.ASprofile(t))/1e-6;


    G = PD.growthrate(S,T,PD.init_dist.boundaries(2:end));


    % Autotimestepsizer based on CFL condition (eq. 24 in Gunawan 2004)
    if flagdt==0

        % Step Sizer only regards the part of the grid that contains
        % (1-mlim) of the distribution (otherwise, for length dependent
        % growth, large L dominate and Dt is small even though there
        % are no crystals at these sizes).
        GI = boundingBoxFinder(fstar(3:end-1),mlim);
        [~,I]   =   max(abs(G(GI(1):GI(2)-1))./Dy(GI(1):GI(2)-1));
        Dt       =   abs(Dy(I)/G(I));


        if length(PD.sol_time)>2
            nexttline   =   PD.sol_time(find(PD.sol_time>t,1,'first'));
            Dttline =   nexttline-t;
        else
            Dttline     =   inf;
        end

        Dt      =   min([Dtlast*1.5 max([min(Dt) PD.sol_time(end)*1e-5]) PD.sol_time(end)-t Dttline 1/10*(PD.sol_time(end)-PD.sol_time(1))]);   % choose minimum of expressions (not too coarse description for plotting purposes)
        Dtlast = Dt;
        fstar0 = fstar;
    end %flagdt

    t           =   t+Dt; % update time step

    % mini failsafe
    if PD.sol_time(end)-t     <   1e-12
        t   =   PD.sol_time(end);
    end

    %% Growth     
    if abs(c-cs)>eps 
        fstar = hiResGrowth(fstar,G,Dt,Dy,GI);
    else
        fstar = fstar0;
    end
    
    %% Nucleation
    if  c>cs
        if nargin(PD.nucleationrate) == 3
            dist = Distribution(y,fstar(3:end-1));
            J = PD.nucleationrate(S,T,dist);
        else
            J = PD.nucleationrate(S,T);
        end
        fstar(3)  = fstar(3) + J/Dy(1)*Dt;
    end

    %% Calculation of concentration loss due to nucleation and growth
    Deltac = sum((fstar(3:end-1)-fstar0(3:end-1)).*y.^3.*...
        Dy(:)*PD.rhoc.*PD.kv);
    
    %% Dilution
    fstar = fstar - Q * fstar/m * Dt; % dilution due to addition of AS
    Deltac = Deltac + Q*c/m*Dt;
    
    
    c_dummy    =    c-Deltac;
    T_dummy    =    PD.Tprofile(t);
    cs_dummy   =    PD.solubility(T_dummy);     %Solubility

    
    
%% finishing        

    DeltaS          =   c_dummy./cs_dummy-c./cs;
    crel = abs((c_dummy-c)/c);
    csrel = abs((cs_dummy-cs)/cs);

    if t<=PD.sol_time(end)
        % Check if result is (superficially) reasonable

        if  (c_dummy>0 && sum(-fstar(fstar<0))<sum(fstar(fstar>0))*1e-2 &&...
                abs(DeltaS)<Stol && ((DeltaS>=0 && (crel<ctol && csrel <ctol)) || (DeltaS<0) && (crel<ctol*10 && csrel <ctol)) || ...
                flagdt > 50)

            if Dt<1e4*eps
                fstar(fstar<0)=0; % hackebeil
            end

            c    =   c_dummy;
            T    =  T_dummy;

            TIME    =   [TIME t]; %#ok<*AGROW>
            tcount  =   tcount+1;   
            flagdt  =   0;

            F = arrayCrop(fstar,[3;sfx-1])';
            Y(tcount,:) = [F(:)' c];

        else
            % Use a smaller timestep and repeat everything
%                 keyboard
            fstar   =   fstar0;
            t       =   t-Dt;
            flagdt  =   flagdt+1;
            Dt      =   Dt/3;

        end

    end


end

end



