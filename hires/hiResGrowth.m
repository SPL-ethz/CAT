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

