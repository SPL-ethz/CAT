function [TIME,Y] = hiResSlave(PD,varargin)
%% n-dimensional High Resolution Flux limited PBE Solver
% This algorithm solves the PDE arising from the population balance
% equation in 1, 2 or 3 dimensions. The correct computation method is
% chosen depending on the shape of the input. The input structure 'input'
% contains information on the experimental and numerical setup (including
% kinetic data, initial distribution and concentration, etc.). 
%
% This code was written together with the HiRes_nD mediator function, so
% check there if you are uncertain about the meaning of some variables.
%
% Dave Ochsenbein, 23.12.2011

%% Setup and Preparation
mlim = 0; % anteil von partikeln, welcher fuer die simulation vernachlaessigt werden darf --> HAS TO BE 0 WHEN WORKING WITH HEAVYSIDE FUNCTION IN F0
% keyboard
if ~isfield(input,'offline')
    input = hiResOLPreparator(input);
end

% Die Idee ist spaeter das hier alles abzuschaffen und nur noch input zu
% passen
f               = PD.init_dist.F;
fstar           = f;

c            =   PD.init_conc;  


TIME = PD.sol_time(1);
%% Integration
flagdt=0;
t=input.offline.t0;tcount=1;

% keyboard
sfx = size(input.offline.f);sfx(sfx==1)=[];
F = arrayCrop(input.offline.f,[3;sfx-1]);

    while t<input.exp.ttot
        T = PD.Tprofile(t);
        cs = PD.solubility(T);
        S = c/cs;
        % Find Growth Rates along all dimensions
        if c>cs+eps
            G = PD.growthrate(S,T,PD.init_dist.boundaries);
        end
        
        % Autotimestepsizer based on CFL condition (eq. 24 in Gunawan 2004)
        if flagdt==0
            
            % Step Sizer only regards the part of the grid that contains
            % (1-mlim) of the distribution (otherwise, for length dependent
            % growth, large L dominate and Dt is small even though there
            % are no crystals at these sizes).
            GI = boundingBoxFinder(fstar(3:end-1),ndim,mlim);
            

            [~,I]   =   max(abs(G(GI(1):GI(2)))./Dx);

            Dt       =   abs(Dx/G(I));

            
            if isfield(input.exp,'tline')
                nexttline   =   input.exp.tline(find(input.exp.tline>t,1,'first'));
                Dttline =   nexttline-t;
            else
                Dttline     =   inf;
            end
       
            Dt      =   min([2000 max([min(Dt) input.exp.ttot*1e-5]) input.exp.ttot-t Dttline]);   % choose minimum of expressions (not too coarse description for plotting purposes)
            
%             
        end %flagdt

        t           =   t+Dt; % update time step

        % mini failsafe
        if input.exp.ttot-t     <   1e-12
            t   =   input.exp.ttot;
        end

        Dtstar = Dt;
        
        %% Growth 
        fstar0 = fstar;
        if abs(c-cs)>eps 
            fstar(:,:,:,:,1) = hiResGrowth(fstar(:,:,:,:,1),G,Dtstar,Dx,GI);
            %% Calculation of concentration at next timestep (Component 1, Growth only)

            DeltacGrowth = sum((fstar(3:end-1,:,:,:,1)-fstar0(3:end-1,:,:,:,1))...
                .*Dx...
                *PD.rhoc.*PD.kv);
        else
            fstar = fstar0;
            DeltacGrowth = 0;
        end
        
        c_dummy    =   c-DeltacGrowth;
        T_dummy    = PD.Tprofile(t);
        cs_dummy   =   PD.solubility(T_dummy);     %Solubility
        
        %% Nucleation
        if strcmp(input.setup.J,'on') && c>cs
            f0(3)  = B/Dx;
            
            c_dummy = c_dummy + sum(DeltaCNuc);
        end
        
%% finishing        

        DeltaS          =   abs(c_dummy./cs_dummy-c./cs); 
        DeltaG  = abs(max(PD.growthrate(S,T,PD.init_dist.boundaries))-max(PD.growthrate(S_dummy,T_dummy,PD.init_dist.boundaries)))/max(PD.growthrate(S_dummy,T_dummy,PD.init_dist.boundaries));
        
        if t<=input.exp.ttot
            % Check if result is (superficially) reasonable

            if  (sum(sum(sum(-fstar(fstar<0))))<sum(sum(sum(fstar(fstar>0))))*1e-2 &&...
                    c_dummy>0 && DeltaS<0.1 && DeltaG<0.1|| ...
                    flagdt > 50)
                
                if Dt<1e4*eps
                    fstar(fstar<0)=0; % hackebeil
                end
                
                c    =   c_dummy;
                
                TIME    =   [TIME t]; %#ok<*AGROW>
                tcount  =   tcount+1;   
                flagdt  =   0;
                
                F = arrayCrop(fstar,[3;sfx-1]);
                Y(tcount,1:end-1) = F(:)';
                Y(tcount,end) = c;

            else
                % Use a smaller timestep and repeat everything
                fstar   =   fstar0;
                t       =   t-Dt;
                flagdt  =   flagdt+1;
                Dt      =   Dt/2;
                                                      
            end

        end
       
        
    end


end



