% file f_odesolve.m
% brief contains function f_odesolve

% brief This function computes time marching for the ODE and PDE system of
% equations computed as part of the IMR solver. The function features known
% integration functions: ODE15, ODE23tb (most stable), and ODE45.
function [t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin)
    
    if method == 15
        if divisions == 0
            options = odeset();
        else
            options = odeset('MaxStep',tfin/divisions,'RelTol',1e-8);
        end
        [t,X] = ode15s(bubble,tspan,init,options);
    elseif method == 23
        if divisions == 0
            options = odeset();
        else
            options = odeset('MaxStep',tfin/divisions,'RelTol',1e-8,'AbsTol',1e-8);
        end
        [t,X] = ode23tb(bubble,tspan,init,options);
    elseif method == 45
        if divisions == 0
            options = odeset('NonNegative',1,'AbsTol',1e-8,'RelTol',1e-8);
        else
            options = odeset('NonNegative',1,'MaxStep',tfin/divisions);
        end
        [t,X] = ode45(bubble,tspan,init,options);
    else
        if divisions == 0
            options = odeset('NonNegative',1);
        else
            options = odeset('NonNegative',1,'MaxStep',tfin/divisions,'RelTol',1e-8);
        end
        [t,X] = ode45(bubble,tspan,init,options);
    end
    
end
