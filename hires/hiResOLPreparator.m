function [finput] = hiResOLPreparator(finput)

% Cut of final time from time points
if ~isfield(finput.exp,'ttot') || (isfield(finput.exp,'tline') && ~isscalar(finput.exp.tline))
    finput.exp.ttot  =   finput.exp.tline(end);
    finput.offline.t0              =   finput.exp.tline(1);
else
    finput.offline.t0              =   0;
end

% Information on the finput.offline.PSD
finput.offline.PSD     =   struct('xb',[],'Dx',[],'Dxprod',[],'xp',[],'F',[]);
finput.offline.PSD(1).xb  =   finput.num.xb(1);

% Find out dimensionality
finput.offline.fields  =   fieldnames(finput.offline.PSD.xb);
finput.offline.ndim    =   numel(finput.offline.fields);

% Grid stuff (pivot lengths, Delta x and include ghost points)
finput.offline.PSD(1).Dxprod = 1;
finput.offline.sizeTot = zeros(1,finput.offline.ndim);
for i=1:finput.offline.ndim
    xb_loc              =   finput.offline.PSD(1).xb.(finput.offline.fields{i});
    finput.offline.PSD(1).Dx.(finput.offline.fields{i})  =   xb_loc(2)-xb_loc(1); % no support for geogrid
    finput.offline.PSD(1).xp.(finput.offline.fields{i})  =   (xb_loc(2:end)+xb_loc(1:end-1))/2;
    size_loc            =   size(finput.offline.PSD(1).xp.(finput.offline.fields{i}));
    finput.offline.sizeTot(i)         =   size_loc(2)+3;
    
    finput.offline.PSD(1).Dxprod = finput.offline.PSD(1).Dxprod * finput.offline.PSD(1).Dx.(finput.offline.fields{i}) ;
end

% keyboard
if isscalar(finput.offline.sizeTot)
    finput.offline.sizeTot = [finput.offline.sizeTot 1];
end

% Copy f0 into (finput.offline.ndim+1)-dimensional array finput.offline.f (extension into time domain)
finput.offline.f = finput.exp.f0;
try
    finput.offline.f = padarray(finput.offline.f,[ones(1,finput.offline.ndim)*2 zeros(1,5-finput.offline.ndim)]);
    finput.offline.f(end,:,:,:,:) = [];
    if finput.offline.ndim>1
        finput.offline.f(:,end,:,:,:) = [];
    end
    if finput.offline.ndim>2
        finput.offline.f(:,:,end,:,:) = [];
    end
catch  %#ok<CTCH>
    finput.offline.f                       =   zeros([finput.offline.sizeTot 1]);  % extend into time domain
    if finput.offline.ndim==1
        finput.offline.f(3:end-1,1,1,1,:)                =   finput.exp.f0;
    elseif finput.offline.ndim==2
        finput.offline.f(3:end-1,3:end-1,1,1,:)          =   finput.exp.f0;
    elseif finput.offline.ndim==3
        finput.offline.f(3:end-1,3:end-1,3:end-1,1,:)    =   finput.exp.f0;
    end
end


T(1)            =   finput.exp.PWCT(1,1); % Initial Temperature

c(1)            =   finput.exp.c0;  
if strcmp(finput.setup.transformation,'on')
    cs(:,1)           =   [finput.data.cs{1}(T(1));...
                            finput.data.cs{2}(T(1))]; % Solubility
else
    cs(:,1)           =   finput.data.cs{1}(T(1));
end

if strcmp(finput.setup.transformation,'on')
    finput.offline.PSD(2).xb  =   finput.num.xb(2);


    for i=1:finput.offline.ndim
        xb_loc              =   finput.offline.PSD(2).xb.(finput.offline.fields{i});
        finput.offline.PSD(2).Dx.(finput.offline.fields{i})  =   xb_loc(2)-xb_loc(1); % no support for geogrid
        finput.offline.PSD(2).xp.(finput.offline.fields{i})  =   (xb_loc(2:end)+xb_loc(1:end-1))/2;
%         size_loc            =   size(finput.offline.PSD(2).xp.(finput.offline.fields{i}));
%         finput.offline.sizeTot(i)         =   size_loc(2)+3;
    end

end
xpArr = xpArrayGenerator(finput.offline.PSD(1).xp,finput.num.ngrid);
finput.offline.xpArr = xpArr;
    
if finput.offline.ndim == 1
    finput.offline.V = finput.data.kv*xpArr.dim1.^3;
elseif finput.offline.ndim == 2
    finput.offline.V = finput.data.kv*xpArr.dim1.*xpArr.dim2.^2;
elseif finput.offline.ndim == 3
    finput.offline.V = finput.data.kv*xpArr.dim1.*xpArr.dim2.*xpArr.dim3;
end
%% Agglo or breakage Preparation
if strcmp(finput.setup.beta,'on') || (strcmp(finput.setup.kbr,'on') && strcmp(finput.setup.millMode,'off'))
    [finput.offline.betaCol,finput.offline.laes,finput.offline.V,finput.offline.vtot,finput.offline.Is] = betaColNSA(finput,xpArr);
else
    finput.offline.betaCol = [];
    finput.offline.laes = [];
    
    finput.offline.vtot = [];
    finput.offline.Is = [];
end

if (strcmp(finput.setup.kbr,'on') || strcmp(finput.setup.kbr,'simple')) && strcmp(finput.setup.millMode,'off')
    [finput.offline.S,finput.offline.BH] = breakageSelection(finput);
elseif strcmp(finput.setup.kbr,'on') && strcmp(finput.setup.millMode,'on')
    [finput.offline.BH1] = breakageSelectionMilling(finput,1);
    [finput.offline.BH2] = breakageSelectionMilling(finput,2);
end

if strcmp(finput.setup.beta,'on')
    
    [finput.offline.tctr0] = tctr0NSA(finput,finput.offline.laes,finput.offline.laes);
     
    finput.offline.fagglo = zeros(finput.offline.sizeTot-3);
    if strcmp(finput.setup.APSD,'on')
                                                                                        
        finput.offline.sourcehandler = @(f0,beta2,Dx,t,Dt,aggloinput,c0,T,finput,G,Lfin) aggloNSA(f0,beta2,Dx,t,Dt,aggloinput,c0,T,finput,G,Lfin);
        
    else
        finput.offline.sourcehandler = @(f0,beta2,Dx,t,Dt,aggloinput,c0,T,finput,G,Lfin) aggloNSAnoAPSD(f0,beta2,Dx,Dt);
    end
    finput.offline.timeResultAgglo.N = zeros(19,1);
    finput.offline.timeResultAgglo.xb = linspace(0,max(max(2*finput.offline.V))^(1/3),20)';
    finput.offline.timeResultAgglo.xp = (finput.offline.timeResultAgglo.xb(1:end-1)+finput.offline.timeResultAgglo.xb(2:end))/2;
else
    finput.offline.timeResultAgglo = [];
    finput.offline.sourcehandler = [];
    finput.offline.tctr0 = [];
end


