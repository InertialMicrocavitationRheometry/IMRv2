function [Tout,Rout,Rddot,Pinf] = f_bubblewall_solver(R_eqd,deltap,kappa,f,...
    We,Ca,rho,po,pv,c,tfinal,force,rmodel,emodel,vmaterial)
% FUNCTION THAT SOLVED KELLER MIKSIS FOR MULTIPLE CASES
% LETTY: YOUR CODE IS A BIT BETTER AT MANAGING VARIABLES, USE YOUR CODE AND
% USE THIS CODE AS  REFERENCE ON HOW TO SIMPLIFY AND DEBUG YOUR CODE
    Ro_eq=R_eqd;                % Bubble Initial Radius
    tf=tfinal;                  % Final time
    c2=c^2;                     % speed of sound squared
    S =(rho*c2*R_eqd)/We;       % Weber #
    eta=(rho*c^2)/Ca;           % eta, SEE JFM PAPER
    pGo=po-pv+(2*S)/R_eqd;      % pressure inside the bubble
    xi = 1/c; 
    if strcmp('RP',rmodel) == 1
        xi = 0;
    elseif strcmp('KM',rmodel) == 1
        xi = 1/c;
    end
    r_dt      = 0;%(pv-po)/(rho*c);       % ICs: initial velocity
    r_initial = [Ro_eq r_dt];             % ICs: initial radius 
    t_span    = [0 tf];                   % time
    options   = odeset('MaxStep', tfinal/1E4,'AbsTol',1E-9);
    % Solving the System of ODEs
    [Tout,Rout] = ode45(@(t,r) bubblewall(t,r,eta,kappa,xi,rho,Ro_eq,...
            S,deltap,f,po,pv,pGo,force,emodel,vmaterial),...
            t_span,r_initial,options);
    [Rddot]     = bubblewall_Rddot(Tout,Rout,eta,kappa,xi,rho,Ro_eq,...
        S,deltap,f,po,pv,pGo,force,emodel,vmaterial);
    [Pinf] = pinf(Tout,f,deltap,force,po);
end
 
function [ dR_dt ] = bubblewall(t,r,eta,k,xi,rho,Req,S,deltap,f,...
    po,pv,pGo,force,emodel,vmaterial)
% KELLER MIKSIS EQUATION NON-LINEAR ODE SOLVER
%
    dRdT2 = bubblewall_Rddot(t,r,eta,k,xi,rho,Req,S,deltap,f,...
    po,pv,pGo,force,emodel,vmaterial);
    dR_dt=[ r(2); dRdT2 ];
end
 
function [ dRdT2 ] = bubblewall_Rddot(t,r,eta,k,xi,rho,Ro,S,deltap,f,...
    po,pv,pGo,force,emodel,vmaterial)
% KELLER MIKSIS EQUATION NON-LINEAR ODE SOLVER
%
    if size(r,2) > 1
        r2 = r(:,2); r1 = r(:,1);
    else
        r2=r(2); r1=r(1);
    end
    r2s=r2.^2;
    pinfy=pinf(t,f,deltap,force,po);
    [mu,mu_inf,mu_o,nc,lambda] = f_carreau(vmaterial,-2*r2./r1);
    nu = mu/rho;
    DEN = (1-r2*xi+4*nu*xi./r1).*r1;
    INE = (-3/2)*(1-r2.*xi/3).*r2s;
    PINF = (1+r2*xi).*((pv-pinfy)/rho);
    PGO = (1+(-3*k+1)*r2*xi).*(pGo/rho).*((Ro./r1).^(3*k));
    if (strcmp('mu_inf',vmaterial) == 1 || strcmp('mu_0',vmaterial) == 1)
        VIS = -4*nu.*r2./r1;
    else
        VIS = -vis(mu_inf,mu_o,nc,lambda,r2,r1)/rho;
    end
%     VIS = 0*VIS;
    SUR = -2*S./(rho.*r1);
    DPA=-((r1/rho)*xi).*pa(t,f,deltap,force);
    % GO THROUGH THE MATH AND SHOW ME THAT YOU CAME UP WITH THESE RESULTS  
    if strcmp('NeoH',emodel)==1
        ELS=-(1./rho).*(1+r2.*xi).*(eta/2).*(5-4*(Ro./r1)-(Ro./r1).^4);
        DELS=-((r1./rho).*xi).*(eta/2).*(4*(Ro./r1.^2).*r2+4*((Ro.^4)./(r1.^5)).*r2);
    elseif strcmp('YangChurch',emodel)==1
        ELS=-(1/rho)*(1+r2*xi).*((4*eta)/3).*(1-(Ro./r1).^3);
        DELS=-((r1./rho).*xi).*(4*eta/3).*r2.*(3*(Ro.^3)./(r1.^4));
    elseif strcmp('LinElas',emodel)==1
        ELS=-(1/rho)*2*eta*(1-(Ro./r1).^2);
        DELS=-((r1/rho)*xi)*4*eta*((Ro.^2)./(r1.^3)).*r2;
    else
        ELS=0; 
        DELS=0;
    end
    dRdT2 = (INE+PINF+PGO+VIS+SUR+ELS+DELS+DPA)./DEN;
end

function dpa = pa(t,f,deltap,force)
% PRESSURE DIFFERENCE CALCULATION
 if strcmp('sine',force)==1
     dpa = deltap*2*pi*f*cos(2*pi*f*t);
 elseif strcmp('mono',force)==1
     dpa = 0;
 end 
end
 
function p = pinf(t,f,deltap,force,po)
% FORCING PRESSURE CALCULATION
 if strcmp('sine',force)==1
     p = deltap*sin(2*pi*f*t);
 elseif strcmp('mono',force)==1
     p = deltap;
 end     
 p = p + po;
end

function S = vis(mu_inf,mu_o,nc,lambda,Rdot,R)
    S = zeros(size(R));
    for i = 1:length(R)
        S(i) = 4*mu_inf*Rdot(i)./R(i)+...
            integral(@(r) carreautau(r,Rdot(i),R(i),mu_inf,mu_o,nc,lambda),R(i),Inf);
    end
end

function tau = carreautau(r,Rdot,R,mu_inf,mu_o,nc,lambda)
    tau = (12*(Rdot.*R.^2)./r.^4).*(mu_o-mu_inf).*(1+...
        lambda.^2.*(4.*Rdot.^2.*R.^4./r.^6)).^((nc-1)/2);    
end