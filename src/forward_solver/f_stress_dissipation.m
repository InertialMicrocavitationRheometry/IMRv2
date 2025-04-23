% file f_stress_dissipation.m
% brief contains function f_stress_dissipation

% brief This function computes the stress dissipation term, \tau : \nabla u.
% The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [taugradu] = f_stress_dissipation(stress,spectral,Req,R,Rdot,Ca,Br, ...
    Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1,ivisco2,fnu,DRe)

% current ammplification factor
Rst = Req/R;
% incompressible condition
x2 = (yT3-1+Rst.^3).^(2/3);
% inversion square of the reference coordinate
ix2 = x2.^-1;
% fourth power of the reference coordinate
x4 = x2.^2;

% no stress
if stress == 0
    taugradu = 0;
    
    % Kelvin-Voigt neo-Hookean solid
elseif stress == 1
    taugradu = 12*(Br/(Re8+DRe*fnu))*(Rdot/R)^2*iyT6 + ...
        2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4);
    
    % Kelvin-Voigt quadratic neo-Hookean solid
elseif stress == 2
    taugradu = 12*(Br/(Re8+DRe*fnu))*(Rdot/R)^2*iyT6 + ...
        2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4) .* ...
        (1+alphax*(x4.*iyT4 + 2*yT2.*ix2 - 3));
    
    % SLS with neo-Hookean solid
elseif stress == 3
    taugradu = 12*(Br/(Re8+DRe*fnu))*(Rdot/R)^2*iyT6 + ...
        2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4);
    
    % SLS with quadratic neo-Hookean solid
elseif stress == 4
    taugradu = 12*(Br/(Re8+DRe*fnu))*(Rdot/R)^2*iyT6 + ...
        2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4) .* ...
        (1+alphax*(x4.*iyT4 + 2*yT2.*ix2 - 3));
    
    % SLS with neo-Hookean solid
elseif stress == 5
    taugradu = 12*(Br/(Re8+DRe*fnu))*(Rdot/R)^2*iyT6 + ...
        2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4);
    
end

% spectral solution approach to compute the mechanical dissipation
if spectral == 1
    taugradu = -2*Br*Rdot./(R*yT3).* ...
        (ZZT*(X(ivisco1)-X(ivisco2)));
end

end
