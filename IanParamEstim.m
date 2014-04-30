function IanParamEstim
% To run this function you need to have downloaded CAT and have copied the
% files SeededGrowthDesupersat_x_T_dave and vanillin_solubility_IR_polyfit
% in the main folder

% load Data
load('SeededGrowthDesupersat_x_T_dave')
load('vanillin_solubility_IR_polyFit');

expii = [2 4]; % indices of experiments in growthExp structure that you want to fit
kitty(length(expii)) = CAT; % create a new array of CATs

% Initial Parameter Guess
p0 = [1.4 0.75 30];

% output vector (concatenate the outputs from different experiments)
y = [];
for ii = 1:length(expii)
    kitty(ii).init_dist = Distribution(growthExp(expii(ii)).PSD.initial.sizeGrid,growthExp(expii(ii)).PSD.initial.q0); %initial distribution
    kitty(ii).init_seed = growthExp(expii(ii)).m_seeds; %seed mass
    kitty(ii).init_massmedium = growthExp(expii(ii)).m_solv; %solvent mass
    kitty(ii).rhoc = 1e-12; %crystal density
    kitty(ii).kv = pi/6; % shape factor

    % Time, temperature and concentration AFTER seeds have been added
    tvec = growthExp(expii(ii)).time;
    Tvec = growthExp(expii(ii)).T;
    xvec = growthExp(expii(ii)).x_vanillin;
    cvec = xvec./(1-xvec);
    tvec = tvec-growthExp(expii(ii)).timeSeedsAdded;
    Tvec(tvec<0) = [];cvec(tvec<0) = [];
    tvec(tvec<0) = [];

    % Fit piecewise-linear function with 3 nodes to data
    PWLfit = slmengine(tvec,Tvec,'degree',1,'knots',3,'decreasing','on');

    % set profiles
    kitty(ii).Tprofile = [PWLfit.knots(:)';PWLfit.coef(:)'];

    % set solubility function (coefficients from Ian)
    kitty(ii).solubility = @(T) polyval([5.94552909168843e-09,-1.28774924731669e-07,5.43215146852858e-06,9.27318853682009e-05,0.00356941579893050],T)./(1-polyval([5.94552909168843e-09,-1.28774924731669e-07,5.43215146852858e-06,9.27318853682009e-05,0.00356941579893050],T));
    kitty(ii).sol_method = 'mp'; % moving pivot
    kitty(ii).init_conc = cvec(1);

    kitty(ii).sol_time = tvec;
    
    y = [y(:); cvec(:)];
    
end
keyboard
options = optimset('tolfun',1e-4);

pstar = fmincon(@(p) objFun(p,kitty,y),p0,[],[],[],[],zeros(size(p0)),[],[],options);

keyboard
end


function f = objFun(p,kitty,y)
fprintf('Parameter Values: p1 =  %4.2f, p2 = %4.2f, p3 = %4.2f \n',p)
% reset growth rate in all CAT objects
for ii = 1:length(kitty)
    kitty(ii).growthrate = @(S,T,y) (S>=1)*p(1)*(S-1)^p(2)*exp(-p(3)/T)*ones(size(y));
end
kitty.solve; % solve it

% construct vector of model outputs (concatenated)
yhat = [];
for ii = 1:length(kitty)
   yhat =  [yhat(:); kitty(ii).calc_conc(:)];
end

% objective function is just the euclidian norm of residuals (square of this is least
% squares)
f = norm(y-yhat);

fprintf('Objective Function Value: f = %4.2g \n',f)

% plot on the go
clf % clear figures
plot(y,'ro','linewidth',1.2)
hold on
plot(yhat,'k-','linewidth',1.2)
drawnow % draw plots before you do something else
keyboard
end