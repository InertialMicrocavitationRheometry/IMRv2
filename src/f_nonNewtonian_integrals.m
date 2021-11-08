function [f,intf,dintf,ddintf] = f_nonNewtonian_integrals(vmodel,U,R)
%F_NONNEWTONIAN_INTEGRALS Summary of this function goes here
    % Setting default values
    f = 0;
    intf = 0;
    dintf = 0;
    ddintf = 0;
    a = 0; nc = 0; lambda = 0;
    gammadot_R   = -2*U/R;
    gammadot_num = -2*U*R*R;
    dgammadot    = -4*U*U*R;
    ddgammadot   = -2*R*R; %the Udot that is calculated in m_cavitation
    abstol = 1E-10; reltol = 1E-10;
    % Setting the parameters and calculating integrals
    switch vmodel
        case 'newtonian'
            return;     
        case 'blood_carreau'
            nc = 0.3568; 
            lambda = 3.313; 
            f = sf_carreau(nc,lambda,gammadot_R);
            I1 = integral(@(r) sf_carreau_d(r,nc,lambda,gammadot_num),...
                R,Inf,'RelTol',reltol,'AbsTol',abstol);
            I2 = integral(@(r) sf_carreau_dd(r,nc,lambda),R,Inf,...
                'RelTol',reltol,'AbsTol',abstol);
            intf = gammadot_num*I1;
            dintf = dgammadot*I2;
            ddintf = ddgammadot*I2;
        case 'blood_carreau_merrill'
            nc = 0.3568; 
            lambda = 3.313; 
            f = sf_carreau(a,nc,lambda,gammadot_R);
        case 'blood_carreau_biro'
            nc = 0.3568; 
            lambda = 3.313;             
            f = sf_carreau(nc,lambda,gammadot_R);
        case 'blood_carreau_skalak'
            nc = 0.3568; 
            lambda = 3.313;             
            f = sf_carreau(nc,lambda,gammadot_R);
        case 'blood_carreau_yasuda'
            f = sf_carreau_yasuda(nc,lambda,gammadot_R);
        case 'blood_powell_eyring'
            f = sf_powell_eyring(nc,lambda,gammadot_R);
        case 'blood_modified_powell_eyring'
            f = sf_modified_powell_eyring(nc,lambda,gammadot_R);
        case 'blood_cross'
            f = sf_cross(nc,lambda,gammadot_R);
        case 'blood_simplified_cross'
            f = sf_simplified_cross(lambda,gammadot_R);
        case 'blood_modified_cross'
            f = sf_modified_cross(a,nc,lambda,gammadot_R);
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%NON NEWTONIAN MODEL EVALUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = sf_powell_eyring(nc,lambda,gammadot_R)
    f = sinh(lambda*gammadot_R)/(lambda*gammadot_R)^nc;
end

function f = sf_modified_powell_eyring(nc,lambda,gammadot_R)
    f = log(lambda*gammadot_R+1)/(lambda*gammadot_R)^nc;
end

function f = sf_cross(a,lambda,gammadot_R)
    f = 1/(1+(lambda*gammadot_R).^a);
end

function f = sf_simplified_cross(lambda,gammadot_R)
    f = 1/(1+(lambda*gammadot_R));
end

function f = sf_modified_cross(a,nc,lambda,gammadot_R)
    f = 1/((1+(lambda*gammadot_R)^nc).^a);
end

function f = sf_carreau(nc,lambda,gammadot_R)
    f = (1+(lambda).^2.*(gammadot_R).^2).^((nc-1)./2);
end

function intf = sf_carreau_d(r,nc,lambda,gammadot_num)
    gammadot = gammadot_num/r^3;
    % additional r in the denominator is from the stress integral
    intf = (gammadot/r)*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
end

function dintf = sf_carreau_dd(r,nc,lambda,gammadot_num)
    % numerator is added later
    gammadot = gammadot_num/r^3;
    % term inside the power, reduces computation
    pref = (1+(lambda).^2.*(gammadot).^2);
    % regular f calculation
    f = pref^((nc-1)./2);
    % derivative of f with respect to gamma only
    dfdgamma = ((nc-1)./2)*pref^((nc-3)./2)*2*lambda^2*gammadot;
    % collecting terms for integration
    dintf = (1/r^4)*(f+gammadot*dfdgamma);
end

function f = sf_carreau_yasuda(a,nc,lambda,gammadot_R)
    f = (1+(lambda).^a.*(gammadot_R).^a).^((nc-1)./a);
end