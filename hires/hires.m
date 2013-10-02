function [TIME,Y] = hiResSlave(PD,varargin)

%% Setup and Preparation
mlim = 0; % anteil von partikeln, welcher fuer die simulation vernachlaessigt werden darf --> HAS TO BE 0 WHEN WORKING WITH HEAVYSIDE FUNCTION IN F0
c            =   PD.init_conc;  


TIME = PD.sol_time(1);
%% Integration
flagdt=0;
t=PD.sol_time(1);tcount=1;

fstar = PD.init_dist.F(:);
x = PD.init_dist.y(:);
fstar = [0;0;fstar;0];
sfx = length(fstar);
F = arrayCrop(fstar,[3;sfx-1]);
Y(1,:) = [F(:)' c];
Dx = diff(PD.init_dist.boundaries);

    while t<PD.sol_time(end)
        % Current solvent + antisolvent mass
        m = PD.init_massmedium+(PD.ASprofile(t)-PD.ASprofile(0));

        % Current mass fraction antisolvent
        xm = PD.ASprofile(t)/m;

        % Current temperature
        T = PD.Tprofile(t);
        
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
            GI = boundingBoxFinder(fstar(3:end-1),1,mlim);
            
            try
            [~,I]   =   max(abs(G(GI(1):GI(2)-1))./Dx(GI(1):GI(2)-1));
            catch
                keyboard
            end
            Dt       =   abs(Dx(I)/G(I));

            
            if length(PD.sol_time)>2
                nexttline   =   PD.sol_time(find(PD.sol_time>t,1,'first'));
                Dttline =   nexttline-t;
            else
                Dttline     =   inf;
            end
       
            Dt      =   min([2000 max([min(Dt) PD.sol_time(end)*1e-5]) PD.sol_time(end)-t Dttline 1/10*(PD.sol_time(end)-PD.sol_time(1))]);   % choose minimum of expressions (not too coarse description for plotting purposes)
            
%             
        end %flagdt

        t           =   t+Dt; % update time step

        % mini failsafe
        if PD.sol_time(end)-t     <   1e-12
            t   =   PD.sol_time(end);
        end

        Dtstar = Dt;
        
        %% Growth 
        fstar0 = fstar;
        if abs(c-cs)>eps 
            fstar = hiResGrowth(fstar,G,Dtstar,Dx,GI);
            
            %% Calculation of concentration loss due to growth (Growth only)
            DeltacGrowth = sum((fstar(3:end-1)-fstar0(3:end-1)).*x.^3.*...
                Dx(:)*PD.rhoc.*PD.kv);
            
            fstar = fstar - Q * fstar/m * Dt; % dilution due to addition of AS
            DeltacGrowth = DeltacGrowth + Q*c/m*Dt;
        else
            fstar = fstar0;
            DeltacGrowth = 0;
        end
        
        c_dummy    =    c-DeltacGrowth;
        T_dummy    =    PD.Tprofile(t);
        cs_dummy   =    PD.solubility(T_dummy);     %Solubility
        
        %% Nucleation
        if  c>cs
            J = PD.nucleationrate(S,T,fstar(3:end-1));
            fstar(3)  = fstar(3) + J/Dx(1)*Dt;
            DeltaCNuc = J*x(1)^3*PD.kv*PD.rhoc*Dt;
            
            c_dummy = c_dummy - sum(DeltaCNuc);
        end
        
%% finishing        

        DeltaS          =   c_dummy./cs_dummy-c./cs;
        S_dummy = c_dummy./cs_dummy;
        if max(PD.growthrate(S,T,PD.init_dist.boundaries)) > eps && max(PD.growthrate(S_dummy,T_dummy,PD.init_dist.boundaries))
            DeltaG  = abs(max(PD.growthrate(S,T,PD.init_dist.boundaries))-max(PD.growthrate(S_dummy,T_dummy,PD.init_dist.boundaries)))/max(PD.growthrate(S_dummy,T_dummy,PD.init_dist.boundaries));
        else
            DeltaG = 0;
        end
        
        if t<=PD.sol_time(end)
            % Check if result is (superficially) reasonable
            
            if  (sum(-fstar(fstar<0))<sum(fstar(fstar>0))*1e-2 &&...
                    c_dummy>0 && DeltaS<0.05 && DeltaG<0.05|| ...
                    flagdt > 50)
                
                if Dt<1e4*eps
                    fstar(fstar<0)=0; % hackebeil
                end
                
                c    =   c_dummy;
                
                TIME    =   [TIME t]; %#ok<*AGROW>
                tcount  =   tcount+1;   
                flagdt  =   0;
                
                F = arrayCrop(fstar,[3;sfx-1])';
                Y(tcount,:) = [F(:)' c];

            else
                % Use a smaller timestep and repeat everything
                fstar   =   fstar0;
                t       =   t-Dt;
                flagdt  =   flagdt+1;
                Dt      =   Dt/2;
%                 keyboard
                                                      
            end

        end
       
        
    end

if length(PD.sol_time)>2
    [~,I] = intersect(TIME,PD.sol_time);
    TIME = TIME(I);
    Y = Y(I,:);
end
end



