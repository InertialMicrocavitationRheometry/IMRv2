% file f_initial_stress_calc.m
% brief contains function f_initial_stress_calc

% brief This function computes the stress value at the bubble maximum
% radius by using shooting method to determine the initial stress. The
% model does not contain heat or mass transfer.
function [S0] = f_init_stress(Req, Re, Ca, De, We, CL, Pv_star)
    sopts = optimset('TolFun', 1e-8);
    Rdotic = fminsearch(@(Rdoti) abs(1-f_maxwell_growth_iter(Rdoti, Req, Re, Ca, De, We, CL)),1, sopts);
    tf_nd = 2;
    tspan = [0, tf_nd];
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    z0 = [Req, Rdotic, 0];
    [~,z] = ode23tb(@(t,z) f_RP_bub_collapse_comp(z, Req, Re, Ca, De, We, CL),tspan,z0, opts);
    [~,Imax] = max(abs(z(:,1)));
    S0 = z(Imax,3);
    
    function Rmax = f_maxwell_growth_iter(Rdoti, Req, Re, Ca, De, We, CL)
        z0 = [Req, Rdoti, 0];
        % final time
        tf_nd = 2;
        tspan = [0, tf_nd];
        opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
        [~,z] = ode23tb(@(t,z) f_RP_bub_collapse_comp(z, Req, Re, Ca, De, We, CL),tspan,z0, opts);
        [Rmax, ~] = max(abs(z(:,1)));
    end
    
    % non-dimensional rayleigh plesset solver
    function [Xdot] = f_RP_bub_collapse_comp(z, Req, Re, Ca, De, We, CL)
        
        % bubble radius R(t)
        R = X(1);
        % bubble wall velocity R'(t)
        U = z(2);
        % stress integral
        S = z(3);
        
        Pgo = 1-Pv_star+1/(We*Req);
        Pb = Pv_star+Pgo*(Req/R)^(3);
        Pbdot = -3*Pgo*(Req/R)^(3)*U/R;
        if De == 0 || De == Inf
            S = (-0.5/Ca*(5-4*Req/R - (Req/R)^4) - 4*U/(Re*R))*R^3;
            S_prime = (2*U*((Req/R)^2 - (Req/R)^5)/Ca - 4/Re*U^2/R^2)*R^3;
            U_prime = ((1+U/CL)*(Pb - 1 -1/(We*R) + S/z1^3) + R/CL*(Pbdot + 1/We*U/R^2+ ...
                (S_prime-3*U*R^2*S/R^3)/R^3) ...
                - 1.5*(1-U/(3*CL))*U^2)/((1-U/CL)*R + 4/Re/CL);
        else % SLS or Maxwell
            S_prime = 1/De*(-4*U/(Re*R)-2*De/Ca*U/R*((Req/R)^4+Req/R)- ...
                (S+1/2*1/Ca*(5-(Req/R)^4-4*Req/R)));
            U_prime = ((1+U/CL)...
                *(pB - 1 + S -1/(We*R)) + R/CL*(Pbdot + 1/We*U/R^2 + S_prime) ...
                - 1.5*(1-U/(3*CL))*U^2)/((1-U/CL)*R + 4/Re/CL);
        end
        Xdot = [U;
        U_prime;
        S_prime];
    end
    
end
