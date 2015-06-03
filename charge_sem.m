%%% calculate charge vector
function [y]=charge(Efmn, Efmp, Egh, vF, dsp_int, flag_prime)
%%% physical constant
global q hbar kBT N3D
n0vec = q^2/pi/hbar^2./vF.^2;

if kBT==0 % T=0, degenerate carrier statistics
    if flag_prime==0    % calculate charge density, apply to E-k of E=sqrt((hbar*vF*k)^2-(Eg)^2) at T=0
        y=q*n0vec.*sign(Efm).*(Efm.^2-Egh.^2).*(abs(Efm)>Egh);
        
    elseif flag_prime==1    % calculate derivitive with regard to En (or -Efm)
        y=-2*q*n0vec.*sign(Efm).*Efm.*(abs(Efm)>Egh);
        
    end
else
    zetan = ( Efmn - Egh)./kBT;     % zetan = Ef - Ec = Ef - Em - Eg/2
    zetap = (-Efmp - Egh)./kBT;     % zetap = Ev - Ef = Em - Ef - Eg/2
    
    if flag_prime==0        
        % calculate charge density, apply to E-k of E=sqrt((hbar*vF*k)^2-(Eg)^2) at T=0
        y = 2*q*n0vec.*(kBT^2*FDjx(zetan,1) - kBT^2*FDjx(zetap,1));
        %%% 2 for valley
    elseif flag_prime==1     
        % calculate derivitive with regard to En (or -Efm)
        y = -2*q*n0vec*(1/kBT).*(kBT^2*FDjx(zetan,0)+ kBT^2*FDjx(zetap,0));  
        %%%Need to be checked
    end
    
%     if flag_prime==0        
%         y(2:end-1) = q*N2D.*(FDjx(zetan(2:end-1),0) - FDjx(zetap(2:end-1),0));
%         %%% no valley degeneracy
%     elseif flag_prime==1     
%         y(2:end-1) = -q*N2D*(1/kBT).*(FDjx(zetan(2:end-1),-1) + FDjx(zetap(2:end-1),-1));  
%         %%%Need to be checked
%     end

    if flag_prime==0        
        y(2:end-1) = q*dsp_int*N3D.*(FDjx(zetan(2:end-1),1/2) - FDjx(zetap(2:end-1),1/2));
        %%% no valley degeneracy
    elseif flag_prime==1     
        y(2:end-1) = -q*dsp_int*N3D*(1/kBT).*(FDjx(zetan(2:end-1),-1/2) + FDjx(zetap(2:end-1),-1/2));  
        %%%Need to be checked
    end
      
end

