function [Udot] = f_radial_eq(radial, p, pVap, pf8, pf8dot, iWe, R, U, J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA )

    % Rayleigh-Plesset 
    if radial == 1
        Udot = (p + pVap - 1 - pf8 - iWe/R + J - 1.5*U^2)/R;

    % Keller-Miksis in pressure        
    elseif radial == 2
        Udot = ((1+U./Cstar)*(p - 1 - pf8 - iWe/R + J) ...
            + R/Cstar*(pdot + iWe*U/R^2 +  JdotX - pf8dot) ...
            - 1.5*(1-U/(3*Cstar))*U^2)/((1-U/Cstar)*R + JdotA/Cstar);     
        % Udot = ((1+U./Cstar)*(p + pVap - 1 - pf8 - iWe./R + J) ... % check to see if pVap should be included 
        %     + R./Cstar.*(pdot + iWe.*U./R.^2 +  JdotX - pf8dot) ...
        %     - 1.5.*(1-U./(3.*Cstar)).*U.^2)./((1-U./Cstar).*R + JdotA./Cstar);      

    % Keller-Miksis in enthalpy
    elseif radial == 3    
        hB = (sam/no)*(((p - iWe/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p - iWe/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/Cstar)*(hB - pf8) - R/Cstar*pf8dot ...
            + R/Cstar*hH*(pdot + iWe*U/R^2 + JdotX) ...
            - 1.5*(1 - U/(3*Cstar))*U^2)/((1 - U/Cstar)*R + JdotA*hH/Cstar);

    % Gilmore equation
	elseif radial == 4
        hB = sam/no*(((p - iWe/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p - iWe/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/Cstar)*(hB - pf8) - R/Cstar*pf8dot ...
            + R/Cstar*(hB + hH*(pdot + iWe*U/R^2 + JdotX)) ...
            - 1.5*(1 - U/(3*Cstar))*U^2) / ((1 - U/Cstar)*R + JdotA*hH/Cstar);
        
    end   

end