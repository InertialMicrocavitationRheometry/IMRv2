%%%%% FUNCTION :: BUBBLE WALL SOLVER
%%%%% INPUTS   :: Input are read as an array called inputs
% odesolve : Matlab ODE solver to use 
% tfinal   : final time of the simulation
% gov_eq   : governing equation
% Ro       : initial radius
% rho      : liquid density
% kappa    : ratio of specific heats
% p_atm    : atmospheric pressure at STP
% pv       : water vapor pressure at STP
% We       : Weber number
% vmodel   : viscosity constitutive equation
% Re       : Reynolds number
% a        : speed of sound
% pfucn    : forcing function
% pdelta   : forcing pressure amplitude
% pf       : forcing pressure frequency
% emodel   : solid material (if any) constitutive equation
% Ca       : Cauchy number
function [Tout,Rout,Rddot,Pinf,DPinf] = f_bubblewall_solver(inputs)
    % reading in the inputs for the calculation
    odesolve = inputs(1);
    tfinal = inputs(2);
    gov_eq = inputs(3);
    Ro     = inputs(4);
    rho    = inputs(5);
    kappa  = inputs(6);
    p_atm  = inputs(7);
    pv     = inputs(8); 
    We     = inputs(9);
    vmodel = inputs(10);
    Re     = inputs(11);
    a      = inputs(12);
    pfunc  = inputs(13);
    pdelta = inputs(14);
    pfreq  = inputs(15);
    emodel = inputs(16);
    Ca     = inputs(17);
    % setting the governing equation, let's set up the KM (w/ enthalpy)
	if strcmp('RP',gov_eq) == 1
        xi = 0;
    elseif strcmp('KM',gov_eq) == 1
        xi = 1/c;
    end
    a2=a^2;                     % speed of sound squared
    S =(rho*c2*R_eqd)/We;       % Weber #
    eta=(rho*c^2)/Ca;           % eta, SEE JFM PAPER
    pGo=po-pv+(2*S)/R_eqd;      % pressure inside the bubble

%     r_dt      = (pv-po)/(rho*c);     % ICs: initial velocity
    r_dt      = 0;                     % ICs: initial velocity
    r_initial = [Ro r_dt];             % ICs: initial radius 
    t_span    = [0 tfinal]; 
    options   = odeset('MaxStep', 5E-10,'AbsTol',1E-12,'RelTol',1E-8);
    % Solving the System of ODEs
    [Tout,Rout] = ode15s(@(t,r) bubblewall(t,r,eta,kappa,xi,rho,Ro_eq,...
            S,deltap,f,po,pv,pGo,force,emodel,vmaterial),...
            t_span,r_initial,options);
    [Rddot]     = bubblewall_Rddot(Tout,Rout,eta,kappa,xi,rho,Ro_eq,...
        S,deltap,f,po,pv,pGo,force,emodel,vmaterial);
    [Pinf]  = pinf(Tout,f,deltap,force,po);
    [DPinf] = dpa(Tout,f,deltap,force);
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
    pinfy=f_p8(t,f,deltap,force,po);
    [mu,mu_inf,mu_o,nc,lambda] = f_carreau(vmaterial,-2*r2./r1);
    nu = mu/rho;
    DEN = (1-r2*xi+4*nu*xi./r1).*r1;
    INE = (-3/2)*(1-r2.*xi/3).*r2s;
    PINF = (1+r2*xi).*((pv-pinfy)/rho);
    PGO = (1+(-3*k+1)*r2*xi).*(pGo/rho).*((Ro./r1).^(3*k));
    if (strcmp('mu_inf',vmaterial) == 1 || strcmp('mu_0',vmaterial) == 1)
        VIS = -4*nu.*r2./r1;
    else
        VIS = -f_vis(mu_inf,mu_o,nc,lambda,r2,r1)/rho;
    end
%     VIS = 0*VIS;
    SUR = -2*S./(rho.*r1);
    DPA=-((r1/rho)*xi).*f_dp8(t,f,deltap,force);

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