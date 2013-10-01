function [Phi] = Phifinder(Theta,PhiStr)
%% Phifinder(Theta,PhiStr)
% Given a theta, find the corresponding flux limiter. The user can choose
% between 5 different functions for Phi (cf. Leveque 2002). Default is
% VanLeer's function for the limiter
% Theta = Ratio of differences between 3 cells
% PhiStr = String determining function that should be used
switch PhiStr
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