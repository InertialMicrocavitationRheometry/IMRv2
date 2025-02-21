% file f_max_pre_stress.m
% brief contains function f_call_params

% brief This function sets up Victor to write this.
function [Smaxpred] = f_max_pre_stress(Ro, al_nd, pwv_nd, We, Re, De, Ca, alpha)
    % Compute trc constant
    trc = sqrt(6 * pi) * gamma(11/6) / (5 * gamma(4/3));
    
    % Compute stress-related terms
    fbarst = pi / (sqrt(6) * We * trc);
    
    % Define constants
    gam = 1.4;
    B = 2.1844;
    
    % Compute barotropic correction term
    fbarbc = -(1 - pwv_nd + 1 / (We * Ro)) * B * Ro^(3 * gam) + pwv_nd;
    
    % Compute compressibility correction
    Mc = 1 / al_nd;
    fbarc = -2 * Mc / (Mc + sqrt(Mc^2 + 4 * trc^2));
    
    % Compute viscous correction
    C = 0.46379 + 0.56391 / Re + 5.74916 / Re^2;
    fbarv = -4 * C^2 / (2 * C^2 + sqrt(4 * C^4 + C^2 * Re^2 * trc^2));
    
    % Compute max stress correction
    fbarmax = fbarv + (De / trc) * (fbarv * exp(-trc / De) - fbarv);
    
    % Compute elastic correction
    fbare = (1 / (60 * Ca * gamma(5/6))) * gamma(1/3) * ...
        ((40 * sqrt(pi) * Ro * (1 - 3 * alpha)) + ...
    (120 * (-1 + 2 * Ro^3) * alpha * gamma(7/6)) / (Ro * gamma(2/3)) + ...
        (-50 + 177 * alpha) * gamma(5/6) / gamma(4/3));
    
    % Compute stress loss correction
    fbarsls = fbarmax - fbare;
    
    % Compute total stress sum
    fsum = fbarbc + fbarst + fbarc + fbarsls;
    
    % Compute time constant tg
    tg = (5 * sqrt(pi) * gamma(5/6) - 6 * Ro^(5/2) * gamma(4/3) * hypergeom([1/2, 5/6], 11/6, Ro^3)) / ...
        (5 * sqrt(6 - 6 * fsum) * gamma(4/3));
    
    % Compute predicted maximum stress
    Smaxpred = fbarv * (1 - exp(-tg / De));
    
end
