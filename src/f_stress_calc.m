% file f_stress_calc.m
% brief contains function f_stress_calc

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [J,JdotX,Z1dot,Z2dot] = ...
        f_stress_calc(stress,X,Req,R,Ca,De,Re8,U,alphax,ivisco,id,LAM,zeNO,cdd)
    % TODO Need to add non-Newtonian behavior to JdotX
    % ((1-U/C_star)*R + ...
        %  4/Re8/C_star - 6*ddintfnu*iDRe/C_star);
    
    Z1dot = zeros(-1,1);
    Z2dot = zeros(-1,1);
    
    % radial stretch
    Rst = Req/R;
    % no stress
    if stress == 0
        J = 0;
        JdotX = 0;
        
        % Kelvin-Voigt with neo-Hookean elasticity
    elseif stress == 1
        J = -(5 - 4*Rst - Rst^4)/(2*Ca) - 4/Re8*U/R;
        JdotX = -2*U/R*(Rst + Rst^4)/Ca + 4/Re8*(U/R)^2;
        
        % quadratic Kelvin-Voigt with neo-Hookean elasticity
    elseif stress == 2
        J = (3*alphax-1)*(5 - Rst^4 - 4*Rst)/(2*Ca) - 4/Re8*U/R + ...
            (2*alphax/Ca)*(27/40 + (1/8)*Rst^8 + (1/5)*Rst^5 + ...
        0.5*Rst^2 - 2*R/Req);
        JdotX = (1/R)*((3*alphax-1)/(2*Ca))*(4*Rst^4*U + 4*Rst*U) + ...
            4*(U^2)/(Re8*R^2) - (2*alphax/Ca)*(2*U/Req + Rst^8*U/R + ...
        Req^5*U/R^6 + Rst^2*U/R);
        
        % linear Maxwell, Jeffreys, Zener -- neo-Hookean
    elseif stress == 3
        % extract stress auxiliary variable
        Z1 = X(ivisco);
        J = Z1/R^3 - 4*LAM/Re8*U/R;
        % elastic shift Ze
        Ze = -0.5 * (R^3 / Ca) * (5 - Rst^4 - 4*Rst);
        % simplified ZdotSqNH equation
        ZdotSqNH = -1.5 * (R^2 * U / Ca) * (5 - Rst^4 - 4*Rst) ...
            - 2 * (R^3 * U / Ca) * (Rst^4 / R + Rst);
        % stress auxiliary variable integral derivative
        Z1dot = -(Z1 - Ze) / De + ZdotSqNH + (3 * U / R) * (Z1 - Ze) ...
            + 4 * (LAM - 1) / (Re8 * De) * R^2 * U;
        % stress integral derivative
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re8*U^2/R^2;
        
        % linear Maxwell, Jeffreys, Zener -- quadratic neo-Hookean
    elseif stress == 4
        % extract stress auxiliary variable
        Z1 = X(ivisco);
        J = Z1/R^3 - 4*LAM/Re8*U/R;
        strainhard = (3*alphax - 1) / (2*Ca);
        % simplified Ze equation with decimal fractions
        Ze = R^3 * (strainhard * (5 - Rst^4 - 4*Rst) + ...
            (2 * alphax / Ca) * (0.675 + 0.125 * Rst^8 + ...
        0.2 * Rst^5 + 0.5 * Rst^2 - 2 / Rst));
        % simplified ZdotSqNH equation with decimal fractions
        ZdotSqNH = (3 * R^2 * U * Ze / R^3) + ...
            R^3 * (strainhard * ((4 * Rst^4 * U / R) + (4 * Rst * U)) - ...
        (2 * alphax / Ca) * ((2 * U / (Rst * R)) + (Rst^8 * U / R) + ...
            (Rst^5 * U / R) + (Rst^2 * U / R)));
        % stress auxiliary derivative
        Z1dot = -(Z1-Ze)/De + ZdotSqNH + (3*U/R)*(Z1-Ze) + 4*(LAM-1)/(Re8*De)*R^2*U ;
        % stress integral derivative
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re8*U^2/R^2;
        
        % upper-convected Maxwell, OldRoyd-B
    elseif stress == 5
        % extract stress sub-integrals
        Z1 = X(ivisco);
        Z2 = X(ivisco2);
        % compute new derivatives
        Z1dot = -(1/De - 2*U/R)*Z1 + 2*(LAM-1)/(Re8*De)*R^2*U;
        Z2dot = -(1/De +   U/R)*Z2 + 2*(LAM-1)/(Re8*De)*R^2*U;
        J = (Z1 + Z2)/R^3 - 4*LAM/Re8*U/R;
        JdotX = (Z1dot+Z2dot)/R^3 - 3*U/R^4*(Z1+Z2) + 4*LAM/Re8*U^2/R^2;
        
    elseif stress == -1
        JdotX = cdd*zeNO;
        
    else
        error('stress setting is not available');
    end
    
end
% % Giesekus, PTT, or forced spectral
% elseif spectral
%     % extract stress spectrum
%     c = X(ic);
% d = X(id);
%     % inverse Chebyshev transforms and derivatives
%     [trr,dtrr,t00,dt00] = stressdiff(c,d);
%     % new spectral coefficient derivatives
%     exptau = exp(ptt*Re8*De*(trr + 2*t00));
%     Z1dot = stresssolve(-(exptau/De + zeNO*4*U./(yV.^3*R) ...
    %         + eps3*Re8*trr).*trr ...
%         + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dtrr ...
    %         - 4./yV.^3*((1-(Req/R)^3)/(3*Ca) + U/(Re8*R) ...
%         + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*4*LAM/Re8*(U/R)^2./yV.^6);
%     Z2dot = stresssolve(-(exptau/De - zeNO*2*U./(yV.^3*R) ...
    %         + eps3*Re8*t00).*t00 ...
%         + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dt00 ...
    %         + 2./yV.^3*((1-(Req/R)^3)/(3*Ca) + U/(Re8*R) ...
%         + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*10*LAM/Re8*(U/R)^2./yV.^6);
%     % compute stress integral
%     J = 2*sum(cdd.*(c-d));
%     JdotX = 2*sum(cdd.*(Z1dot - Z2dot));

%  functions called by solver
%
% % stress differentiator
% function [trr,dtrr,t00,dt00] = stressdiff(c,d)
%     if Nv < 650
%         trr = sCA*c;
%         dtrr = sCAd*c;
%         t00 = sCA*d;
%         dt00 = sCAd*d;
%     else
%         [trr,dtrr] = fctdShift(c);
%         [t00,dt00] = fctdShift(d);
%     end
% end
%
% % stress solver
% function s = stresssolve(x)
%     if Nv < 650
%         s = sCI*x;
%     else
%         s = fctShift(x);
%     end
% end
%
% % fast Chebyshev transform
% function a = fctShift(v)
%     v = v(:);
%     v = [0;
%     v;
%     flipud(v(1:Nv-1))];
%     a = real(fft(v))/Nv;
%     a = [a(2:Nv);
%     a(Nv+1)/2];
% end
%
% % fast Chebyshev transform and differentiate
% function [v,w] = fctdShift(a)
%     M = Nv + 1;
%     a = a(:)';
%     dd = Nv*[0 a(1:Nv-1) a(Nv)*2 fliplr(a(1:Nv-1))];
%     v = ifft(dd);
%     v = v(2:M)' - sum(a);
%     n2b = (0:M-2).^2.*dd(1:Nv);
%     cc = imag(ifft([0:M-2 0 2-M:-1].*dd));
%     w = zeros(Nv,1);
%     w(1:Nv-1) = csc(pi/Nv*(1:M-2)).*cc(2:Nv);
%     w(Nv) = sum((-1).^(1:Nv).*n2b)/Nv + 0.5*(-1)^M*Nv*dd(M);
% end
