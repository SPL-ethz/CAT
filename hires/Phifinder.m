function [Phi] = Phifinder(Theta,Phistr)

% Create Flux Limiter Anonymous Function
try
    if strcmpi(Phistr,'vanleer')

        Phi     =   (abs(Theta)+Theta)./(1+abs(Theta));

    elseif strcmpi(Phistr,'koren')
        % Check expression here again, might be a little off (different
        % description of theta?)
        for i=1:length(Theta)
            Phi(i)  =   max([0 min([2.*Theta(i) min([1/3+2/3*Theta(i) 2])])]);
        end

    elseif strcmpi(Phistr,'minmod')
        for i=1:length(Theta)
            Phi(i)  =   max([1 Theta(i)]);
        end

    elseif strcmpi(Phistr,'superbee')
        for i=1:length(Theta)
            Phi(i)  =   max([0 min([1 2*Theta(i)]) min([2 Theta(i)])]);
        end

    elseif strcmpi(Phistr,'MC')
        for i=1:length(Theta)
            Phi(i)  =   max([0 min([2 2.*Theta(i) (1+Theta(i))/2])]);
        end
    end
        
catch
    % Default is VanLeer 
    Phi=(abs(Theta)+Theta)./(1+abs(Theta));
    
end

end