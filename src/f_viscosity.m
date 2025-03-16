% file f_viscosity.m
% brief contains function f_viscosity

% brief This function computes the nonNewtonian integrals and accounts for
% Newtonian, Carreau, Carreau-Yasuda, Powell-Eyring, modified
% Powell-Eyring, Cross, simplified Cross, and modified Cross

function [f,intf,dintf,ddintf] = ...
    f_viscosity(vmodel,Rdot,R,a,nc,lambda)

gammadot_R   = -2*Rdot/R;
gammadot_num = -2*Rdot*R*R;
dgammadot    = -4*Rdot*Rdot*R;
ddgammadot   = -2*R*R; %the Rdotdot that is calculated in m_cavitation
abstol = 1E-8;
reltol = 1E-8;
% Setting the parameters and calculating integrals
switch vmodel
    case 'newtonian'
    case 'carreau'
    % calculating the viscosity slope at r = R
    f = sf_carreau(nc,lambda,gammadot_R);
    % calculating the Leibniz integration rule limit, see Appendix
    % of the manuscript
    h = f*gammadot_R*(Rdot/R);
    % calculating the stress integral for a non-Newtonian model,
    % goes directly into the E term in the Keller-Miksis equation
    I1 = integral(@(r) sf_carreau_d(r,nc,lambda,gammadot_num),...
        R,Inf,'RelTol',reltol,'AbsTol',abstol);
    % calculating the time derivative of the stress integral for
    % a non-Newtonian model, this integral has to two coefficients.
    % One of the terms is in the E_primber term in m_cavitation,
    % the other goes in the denominator of the Rdot_dot solution to
    % account for the Rddot term
    I2 = integral(@(r) sf_carreau_dd(r,nc,lambda,gammadot_num),R,Inf,...
        'RelTol',reltol,'AbsTol',abstol);
    intf = gammadot_num*I1;
    % note the additional h term is to account for the Leibniz
    % integration rule correction
    dintf = dgammadot*I2 - h;
    % second term that is used for the denominator of the Rdot_dot
    % calculation
    ddintf = ddgammadot*I2;
    case 'carreau_yasuda'
    f = sf_carreau_yasuda(nc,lambda,gammadot_R);
    case 'powell_eyring'
    f = sf_powell_eyring(nc,lambda,gammadot_R);
    case 'modified_powell_eyring'
    f = sf_modified_powell_eyring(nc,lambda,gammadot_R);
    case 'cross'
    f = sf_cross(nc,lambda,gammadot_R);
    case 'simplified_cross'
    f = sf_simplified_cross(lambda,gammadot_R);
    case 'modified_cross'
    f = sf_modified_cross(a,nc,lambda,gammadot_R);
end

end

%non-Newtonian model evaluations
function intf = sf_carreau_d(r,nc,lambda,gammadot_num)
    gammadot = gammadot_num./r.^3;
    % additional r in the denominator is from the stress integral
    intf = (gammadot./r).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
end

function dintf = sf_carreau_dd(r,nc,lambda,gammadot_num)
    % numerator is multiplied later
    gammadot = gammadot_num./r.^3;
    % term inside the power, reduces computation
    pre_f = (1+(lambda).^2.*(gammadot).^2);
    % regular f calculation
    f = pre_f.^((nc-1)./2);
    % derivative of f with respect to gammadot only
    dfdgamma = ((nc-1)./2)*pre_f.^((nc-3)./2)*2*lambda^2.*gammadot;
    % collecting terms for integration
    dintf = (1./r.^4).*(f+gammadot.*dfdgamma);
end

% function f = sf_carreau(nc,lambda,gammadot_R)
%     f = (1+(lambda).^2.*(gammadot_R).^2).^((nc-1)./2);
% end

% function f = sf_carreau_yasuda(a,nc,lambda,gammadot_R)
%     f = (1+(lambda).^a.*(gammadot_R).^a).^((nc-1)./a);
% end
%
% function f = sf_cross(a,lambda,gammadot_R)
%     f = 1/(1+(lambda*gammadot_R).^a);
% end
%
% function f = sf_modified_cross(a,nc,lambda,gammadot_R)
%     f = 1/((1+(lambda*gammadot_R)^nc).^a);
% end
%
% function f = sf_modified_powell_eyring(nc,lambda,gammadot_R)
%     f = log(lambda*gammadot_R+1)/(lambda*gammadot_R)^nc;
% end
%
% function f = sf_powell_eyring(nc,lambda,gammadot_R)
%     f = sinh(lambda*gammadot_R)/(lambda*gammadot_R)^nc;
% end
%
% function f = sf_simplified_cross(lambda,gammadot_R)
%     f = 1/(1+(lambda*gammadot_R));
% end
