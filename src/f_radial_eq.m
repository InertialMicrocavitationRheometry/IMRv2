% file f_radial_eq.m
% brief contains function f_radial_eq

% brief This function features the spherical radial bubble dynamics
% equations used to compute the bubble wall acceleration. The radial models
% Rayleigh-Plesset, Keller-Miksis in pressure and enthalpy, and Gilmore are
% available in the solver.
function [Rddot] = f_radial_eq(radial, p, pdot, pVap, pf8, pf8dot, iWe, ...
        R, Rdot, J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA )
    
    % Rayleigh-Plesset
    if radial == 1
        Rddot = (p + pVap - 1 - pf8 - iWe/R + J - 1.5*Rdot^2)/R;
        
        % Keller-Miksis in pressure
    elseif radial == 2
        Rddot = ((1+Rdot./Cstar)*(p + pVap - 1 - pf8 - iWe/R + J) ...
            + R/Cstar*(pdot + iWe*Rdot/R^2 +  JdotX - pf8dot) ...
        - 1.5*(1-Rdot/(3*Cstar))*Rdot^2)/((1-Rdot/Cstar)*R + JdotA/Cstar);
        
        % Rddot = ((1+Rdot./Cstar)*(p + pVap - 1 - pf8 - iWe./R + J) ... % check to see if pVap should be included
        %     + R./Cstar.*(pdot + iWe.*Rdot./R.^2 +  JdotX - pf8dot) ...
            %     - 1.5.*(1-Rdot./(3.*Cstar)).*Rdot.^2)./((1-Rdot./Cstar).*R + JdotA./Cstar);
        
        % if fdkv == 1
        %     RHS = (1+Rdot/Cstar)...
            %         *(P  + abs(1-1)*Pv -1/(We*R) + J - 1 - Pext)  ...
        %         + R/Cstar*(Pdot+ Rdot/(We*R^2) + Jdot -P_ext_prime );
        %     LHS = (3/2)*(1-Rdot/(3*Cstar))*Rdot^2;
        %     denom = (1-Rdot/Cstar)*R - (R/Cstar)*Sdd;
        %
        %     Rddot = (RHS - LHS)/denom;
        %
        
        % Keller-Miksis in enthalpy
    elseif radial == 3
        hB = (sam/no)*(((p - iWe/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p - iWe/R + GAMa + J))^(1/nstate);
        Rddot = ((1 + Rdot/Cstar)*(hB - pf8) - R/Cstar*pf8dot ...
            + R/Cstar*hH*(pdot + iWe*Rdot/R^2 + JdotX) ...
        - 1.5*(1 - Rdot/(3*Cstar))*Rdot^2)/((1 - Rdot/Cstar)*R + JdotA*hH/Cstar);
        
        % Gilmore equation
    elseif radial == 4
        hB = sam/no*(((p - iWe/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p - iWe/R + GAMa + J))^(1/nstate);
        Rddot = ((1 + Rdot/Cstar)*(hB - pf8) - R/Cstar*pf8dot ...
            + R/Cstar*(hB + hH*(pdot + iWe*Rdot/R^2 + JdotX)) ...
        - 1.5*(1 - Rdot/(3*Cstar))*Rdot^2) / ((1 - Rdot/Cstar)*R + JdotA*hH/Cstar);
        
    end
    
end