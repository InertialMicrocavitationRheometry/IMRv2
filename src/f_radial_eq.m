% file f_radial_eq.m
% brief contains function f_radial_eq

% brief This function features the spherical radial bubble dynamics
% equations used to compute the bubble wall acceleration. The radial models
% Rayleigh-Plesset, Keller-Miksis in pressure and enthalpy, and Gilmore are
% available in the solver.
function [Rddot] = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, ...
    R, Rdot, J, Jdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
    ddintfnu, iDRe )

% Rayleigh-Plesset
if radial == 1
    Rddot = (P - 1 - Pf8 - iWe/R + J - 1.5*Rdot^2)/R;
    
    % Keller-Miksis in pressure
elseif radial == 2
    Rddot = ((1+Rdot./Cstar)*(P - 1 - Pf8 - iWe/R + J) ...
        + R/Cstar*(Pdot + iWe*Rdot/R^2 +  Jdot - Pf8dot) ...
        - 1.5*(1-Rdot/(3*Cstar))*Rdot^2)/((1-Rdot/Cstar)*R + JdotA/Cstar ...
        - 6*ddintfnu*iDRe/Cstar);
    
    % if fdkv == 1
    %     RHS = (1+Rdot/Cstar)...
        %         *(P  + abs(1-1)*Pv -1/(We*R) + J - 1 - Pext)  ...
        %         + R/Cstar*(Pdot+ Rdot/(We*R^2) + Jdot -P_ext_prime );
    %     LHS = (3/2)*(1-Rdot/(3*Cstar))*Rdot^2;
    %     denom = (1-Rdot/Cstar)*R - (R/Cstar)*Sdd;
    %
    %     Rddot = (RHS - LHS)/denom;
    %
    
    % Keller-Miksis in enthalpy with Tait EoS
elseif radial == 3
    hB = (sam/no)*(((P - iWe/R + GAMa + J)/sam)^no - 1);
    hH = (sam/(P - iWe/R + GAMa + J))^(1/nstate);
    Rddot = ((1 + Rdot/Cstar)*(hB - Pf8) - R/Cstar*Pf8dot ...
        + R/Cstar*hH*(Pdot + iWe*Rdot/R^2 + Jdot) ...
        - 1.5*(1 - Rdot/(3*Cstar))*Rdot^2) / ((1 - Rdot/Cstar)*R + ...
        JdotA*hH/Cstar - 6*hH*ddintfnu*iDRe/Cstar);
    
    % Gilmore equation with Tait EoS
elseif radial == 4
    Pb = P - iWe/R + GAMa + J;
    rho = (Pb/sam)^(1/nstate);
    C = sqrt(nstate*Pb/rho);
    hB = sam/no*((Pb/sam)^no - 1);
    hH = (sam/Pb)^(1/nstate);
    Rddot = ((1 + Rdot/C)*(hB - Pf8) - R/C*Pf8dot ...
        + R/C*hH*(Pdot + iWe*Rdot/R^2 + Jdot) ...
        - 1.5*(1 - Rdot/(3*C))*Rdot^2) / ((1 - Rdot/C)*R + ...
        JdotA*hH/C - 6*hH*ddintfnu*iDRe/C);
    
    % Keller-Miksis in enthalpy with Mie-Gruneisen EoS
elseif radial == 5
    % nog = (nstate - 1) / 2;
    [~, hB, hH] = f_mie_gruneisen_eos_scalar(P);
    Rddot = ((1 + Rdot/Cstar)*(hB - Pf8) - R/Cstar*Pf8dot ...
        + R/Cstar*hH*(Pdot + iWe*Rdot/R^2 + Jdot) ...
        - 1.5*(1 - Rdot/(3*Cstar))*Rdot^2) / ((1 - Rdot/Cstar)*R + ...
        JdotA*hH/Cstar - 6*hH*ddintfnu*iDRe/Cstar);
    
    % Gilmore with Mie-Gruneisen EoS
elseif radial == 6
    [C, hB, hH] = f_mie_gruneisen_eos_scalar(P);
    Rddot = ((1 + Rdot/C)*(hB - Pf8) - R/C*Pf8dot ...
        + R/C*hH*(Pdot + iWe*Rdot/R^2 + Jdot) ...
        - 1.5*(1 - Rdot/(3*C))*Rdot^2) / ((1 - Rdot/C)*R + ...
        JdotA*hH/C - 6*hH*ddintfnu*iDRe/C);
    
end

% computes mu, rho, enthalpy h, and hH = 1/rho at a scalar pressure P
% using the Mie–Grüneisen EoS
function [C, hB, hH] = f_mie_gruneisen_eos_scalar(pval)
    
    % normalize pressure
    A = pval / (rho0 * C0^2);
    As = A*hugoniot_s;
    As2 = A*hugoniot_s^2;
    
    % quadratic coefficients for solving mu
    a = As2 - nog;
    b = -2*As - 1;
    % discriminant calculation
    d = b^2 - 4*a*A;
    
    %  density strain
    mu = (-b + sqrt(d)) / (2 * a);
    C = Cstar*sqrt(((1+2*nog*mu)*(1-s*mu)^2+2*s*mu*(1+nog*mu)*(1-s*mu))/(1-s*mu)^4);
    
    %   hH  - 1 / ρ(P) (for dot{h} = hH * dot{P})
    hH = 1/(rho0*(1 + mu));
    
    % compute enthalpy numerically: h(P) = ∫_{Pinf}^{P} (1/ρ(p')) dp'
    integrand = @(pp) 1 / f_mie_rho_from_p_scalar(pp);
    % enthalpy: h(P) = ∫_{Pinf}^P (1/ρ(p')) dp'
    hB = integral(integrand, 1, pval, 'AbsTol', 1e-8, 'RelTol', 1e-8);
    
end

% function to compute rho(p) from Mie–Grüneisen EOS (scalar p)
function rho = f_mie_rho_from_p_scalar(p)
    A = p / Cstar^2;
    a = A*hugoniot_s^2 - nog;
    b = -2*A*hugoniot_s - 1;
    d = b^2 - 4 * a * A;
    rho = rho0 * (1 + (-b + sqrt(d))/(2*a));
end

end
