function[y] = hiResGrowth(fstarfull,G,Dt,Dx,GI)
   

FI = zeros(size(GI));
FI(1) = max([GI(1)-3 1]);
FI(2) = min([GI(2)+3 length(fstarfull)]);
Dx = [eps;eps;Dx(:); eps];Dx = Dx(FI(1):FI(2));


if any(max(G)>eps)

    fstarfull = permute(fstarfull,[circshift(1:1,[0 1-1]) 1+1:3 4]);

    fstar = fstarfull(FI(1):FI(2),:,:,:);

    % Find Theta
    if any(G>=  0)
        deltaf  =   diff(fstar,1,1);
        Theta   =   arrayDivision(deltaf,1,1);
    else
        fstar(1:2) = repmat(fstar(3),[2 1]); % outflow boundary condition p 131 ff leveque
        deltaf  =   diff(fstar,1,1);
        Theta   =   arrayDivision(deltaf,1,2);
        Theta   =   circshift(Theta,[-1 0]);
    end

    Theta(isinf(Theta))     =   0;
    Theta(isnan(Theta))     =   2;

    Phi                     =   Phifinder(Theta,'vanleer');    % Flux limiter
    Phi(isnan(Phi))         = 0;

    % Propagate Distribution based on G setup

    if length(unique(G))==1
        Gloc = G(1);
        if min(G)>=0

            f_dummy = fstar(3:end-1)-Dt.*Gloc./Dx(3:end-1).*(fstar(3:end-1)-fstar(2:end-2,:,:,:))-Dt.*Gloc./(2*Dx(3:end-1)).*(1-Dt.*Gloc./Dx(3:end-1)).*...
                ((fstar(4:end,:,:,:)-fstar(3:end-1)).*Phi(2:end,:,:,:)-(fstar(3:end-1)-fstar(2:end-2,:,:,:)).*Phi(1:end-1,:,:,:));

        else
            f_dummy = fstar(3:end-1)-Dt*Gloc./Dx.*(fstar(4:end,:,:,:)-fstar(3:end-1))-Dt*Gloc./(Dx*2).*(1+Dt*Gloc./Dx).*...
                ((fstar(3:end-1)-fstar(2:end-2,:,:,:)).*Phi(1:end-1,:,:,:)-(fstar(4:end,:,:,:)-fstar(3:end-1)).*Phi(2:end,:,:,:));     
        end

    elseif length(unique(G))>1

        Gloc = [0;0;G(:); 0];Gloc = Gloc(FI(1):FI(2));Gloc=Gloc(:);

        if min(Gloc(:))>=0
%             keyboard
            f_dummy=fstar(3:end-1)-Dt./Dx(3:end-1).*(Gloc(3:end-1).*fstar(3:end-1)-Gloc(2:end-2).*fstar(2:end-2,:,:,:))-...           % f(t=t+Dt) saved in a dummy variable for convenience
                            (Dt./(2*Dx(3:end-1)).*Gloc(3:end-1).*(1-(Dt./Dx(3:end-1).*Gloc(3:end-1))).*(fstar(4:end,:,:,:)-fstar(3:end-1)).*Phi(2:end,:,:,:)-...
                            Dt./(2*Dx(2:end-2)).*Gloc(2:end-2).*(1-Dt./Dx(2:end-2).*Gloc(2:end-2)).*(fstar(3:end-1)-fstar(2:end-2,:,:,:)).*Phi(1:end-1,:,:,:));

        else
            f_dummy=fstar(3:end-1,:,:)-Dt./Dx.*(Gloc(4:end,:,:).*fstar(4:end,:,:)-Gloc(3:end-1,:,:).*fstar(3:end-1,:,:))+...           % f(t=t+Dt) saved in a dummy variable for convenience
                        (Dt./(2*PSD.Dx).*Gloc(3:end-1,:,:).*(1+(Dt./Dx.*Gloc(3:end-1,:,:))).*(fstar(4:end,:,:)-fstar(3:end-1,:,:)).*Phi(2:end,:,:)-...
                        Dt./(2*PSD.Dx).*Gloc(2:end-2,:,:).*(1+Dt./Dx.*Gloc(2:end-2,:,:)).*(fstar(3:end-1,:,:)-fstar(2:end-2,:,:)).*Phi(1:end-1,:,:));
%                     
        end

    end

    fstar(3:end-1,:,:)    =   f_dummy;

    fstarfull(FI(1):FI(2),:,:,:) = fstar;

end



y = fstarfull;

