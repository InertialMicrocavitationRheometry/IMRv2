% file f_stress_calc.m
% brief contains function f_stress_calc

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [S,Sdot,Z1dot,Z2dot] = f_stress_calc(stress,X,Req,R,Ca,De,Re8,...
    Rdot,alphax,ivisco1,ivisco2,LAM,zeNO,cdd,intfnu,dintfnu,iDRe)

Z1dot = [];
Z2dot = [];

% radial stretch
Rst = Req/R;

% no stress
if stress == 0
    S = 0;
    Sdot = 0;
    
    % Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 1
    S = -(5 - 4*Rst - Rst^4)/(2*Ca) - 4/Re8*Rdot/R - 6*intfnu*iDRe;
    Sdot = -2*Rdot/R*(Rst + Rst^4)/Ca + 4/Re8*(Rdot/R)^2 - 6*dintfnu*iDRe;
    
    % quadratic Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 2
    S = (3*alphax-1)*(5 - Rst^4 - 4*Rst)/(2*Ca) - 4/Re8*Rdot/R - 6*intfnu*iDRe + ...
        (2*alphax/Ca)*(27/40 + (1/8)*Rst^8 + (1/5)*Rst^5 + ...
        Rst^2 - 2/Rst);
    Sdot = (Rdot/R)*((3*alphax - 1)/(2*Ca))*(4*Rst^4+4*Rst) + ...
        4*(Rdot/R)^2/Re8 - 6*dintfnu*iDRe -...
        2*alphax/Ca*Rdot/R*(Rst^8 + Rst^5 + 2*Rst^2 + 2*Rst^(-1));
    
    % linear Maxwell, Jeffreys, Zener -- neo-Hookean
elseif stress == 3
    % extract stress auxiliary variable
    Z1 = X(ivisco1);
    S = Z1/R^3 - 4*LAM/Re8*Rdot/R - 6*LAM*intfnu*iDRe;
    % elastic shift Ze
    Ze = -0.5*(R^3/Ca)*(5 - Rst^4 - 4*Rst);
    % stress auxiliary variable integral derivative
    Z1dot = -(Z1 - Ze)/De + 4*(LAM - 1)/(Re8*De)*R^2*Rdot;
    % stress integral derivative
    Sdot = Z1dot/R^3 - 3*Rdot/R^4*Z1 + 4*LAM/Re8*(Rdot/R)^2;
    
    % linear Maxwell, Jeffreys, Zener -- quadratic neo-Hookean
elseif stress == 4
    % extract stress auxiliary variable
    Z1 = X(ivisco1);
    S = Z1/R^3 - 4*LAM/Re8*Rdot/R;
    strainhard = (3*alphax - 1) / (2*Ca);
    % simplified Ze equation with decimal fractions
    Ze = R^3 * (strainhard * (5 - Rst^4 - 4*Rst) + ...
        (2 * alphax / Ca) * (0.675 + 0.125 * Rst^8 + ...
        0.2 * Rst^5 + Rst^2 - 2 / Rst));
    % stress auxiliary derivative
    Z1dot = -(Z1-Ze)/De + 4*(LAM-1)/(Re8*De)*R^2*Rdot ;
    % stress integral derivative
    Sdot = Z1dot/R^3 - 3*Rdot/R^4*Z1 + 4*LAM/Re8*Rdot^2/R^2;
    
    % upper-convected Maxwell, OldRoyd-B
elseif stress == 5
    % extract stress sub-integrals
    Z1 = X(ivisco1);
    Z2 = X(ivisco2);
    % compute new derivatives
    Z1dot = -(1/De - 2*Rdot/R)*Z1 + 2*(LAM-1)/(Re8*De)*R^2*Rdot;
    Z2dot = -(1/De +   Rdot/R)*Z2 + 2*(LAM-1)/(Re8*De)*R^2*Rdot;
    % stress integral and derivative
    S = (Z1 + Z2)/R^3 - 4*LAM/Re8*Rdot/R;
    Sdot = (Z1dot+Z2dot)/R^3 - 3*Rdot/R^4*(Z1+Z2) + 4*LAM/Re8*Rdot^2/R^2;
    
elseif stress == -1
    % this is to keep the code clean
    % need to fix for spectral
    Sdot = cdd*zeNO;
    
else
    error('stress setting is not available');
end

end
