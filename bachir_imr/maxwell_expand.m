% Preliminary study for scaling of Maxwell element stress at end of initial
% expansion

% Zhiren - 4/7/2021

function [t1,X1] = maxwell_expand(lam0,t0,tau)

X0 = 1; % Here, X is elastic stretch
trange = [-t0,0];
[t, X] = ode23tb(@expansion, trange, X0);

% For convenience, let's output the whole thing. We can post-process
% outside, as needed.
t1 = t;
X1 = X;

% ----------------------
% We need to embed a function to advance ODE, since ODE23TB can only have a
% fixed structure of input

function dxdt = expansion(t,x)
    
    % This is for a specific assumed function of initial expansion:
    % lambda = 1 + (lambda0 - 1)*(0.5)*(1+cos(pi*t/t0))
    % Which allows transition between 0 and 1, with 0 speed at both ends.
    
    LHS = (1-lam0)*sin(pi*t/t0)*(pi/t0)/(2+(lam0-1)*(1+cos(pi*t/t0)));
    RHS = log(x)/tau;
    
    dxdt = (LHS-RHS)*x;
    

end % Embeded expansion
    
end % Main function
