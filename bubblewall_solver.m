function [ODE] = bubblewall_solver(R_eqd,A,w,We,Ca,rho,po,pv,c,tfinal,rmodel,emodel,vmodel,vmaterial)
% FUNCTION THAT SOLVED KELLER MIKSIS FOR MULTIPLE CASES
% LETTY: YOUR CODE IS A BIT BETTER AT MANAGING VARIABLES, USE YOUR CODE AND
% USE THIS CODE AS  REFERENCE ON HOW TO SIMPLIFY AND DEBUG YOUR CODE
    Ro_eq=R_eqd;                %Bubble Initial Radius
    tf=tfinal;                  %Final time
    c2=c^2;                     % speed of sound squared
    k=1;                        % gas parameter, using air value
    S =(rho*c2*R_eqd)/We;       % Weber #
    eta=(rho*c^2)/Ca;           % eta, SEE JFM PAPER
    pGo=po-pv+(2*S)/R_eqd;      % pressure inside the bubble
    xi = 1/c; 
    if strcmp('RP',rmodel) == 1
        xi = 0;
    elseif strcmp('KM',rmodel) == 1
        xi = 1/c;
    end
    r_dt=0;%(pv-po)/(rho*c);       % ICs: initial velocity
    r_initial=[Ro_eq r_dt];     % ICs: initial radius 
    t_span=[0 tf];              % TIME
    options = odeset('MaxStep', tf/1E4,'AbsTol',1E-9);
    % Solving the System of ODEs
    [ODE]=ode23tb(@(t,r) bubblewall(t,eta,r,k,xi,rho,Ro_eq,...
            S,A,w,po,pv,pGo,emodel,vmodel,vmaterial),t_span,r_initial,options);
end
 
function [ dR_dt ] = bubblewall(t,eta,r,k,xi,rho,Req,S,A,w,...
    po,pv,pGo,emodel,vmodel,vmaterial)
% KELLER MIKSIS EQUATION NON-LINEAR ODE SOLVER
%
    Ro=Req;
    r2=r(2);
    r1=r(1);
    r2s=r2^2;
    pinfy=pinf(t,w,A)+po;
    if strcmp('Carreau',vmodel) == 1
         mu = carreau(vmaterial,r2/r1);
    elseif strcmp('powerlaw',vmodel) == 1
         mu = powerlaw(vmaterial,r2/r1);
    end
    DEN = (1-r2*xi+4*mu*xi/(rho*r1))*r1;
    INE = (-3/2)*(1-r2*((1/3)*xi))*r2s;
    PINF = (1+r2*xi)*((pv-pinfy)/rho);
    PGO = (1+(-3*k+1)*r2*xi)*(pGo/rho)*((Ro/r1)^(3*k));
    VIS = -4*mu*r2/(r1*rho);%(-4*nu*r2)/r1;
    SUR =-(2*S)/(rho*r1);
    DPA=-((r1/rho)*xi)*pa(t,w,A);
    % GO THROUGH THE MATH AND SHOW ME THAT YOU CAME UP WITH THESE RESULTS  
    ELS = 0;
    DELS = 0;
    if strcmp('NeoH',emodel)==1
        ELS=-(1/rho)*(1+r2*xi)*(eta/2)*(5-4*(Ro/r1)-(Ro/r1)^4);
        DELS=-((r1/rho)*xi)*(eta/2)*(4*(Ro/r1^2)*r2+4*((Ro^4)/(r1^5))*r2);
    elseif strcmp('YangChurch',emodel)==1
        ELS=-(1/rho)*(1+r2*xi)*((4*eta)/3)*(1-(Ro/r1)^3);
        DELS=-((r1/rho)*xi)*(4*eta/3)*r2*(3*(Ro^3)/(r1^4));
    elseif strcmp('LinElas',emodel)==1
        ELS=-(1/rho)*2*eta*(1-(Ro/r1)^2);
        DELS=-((r1/rho)*xi)*4*eta*((Ro^2)/(r1^3))*r2;
    else
        ELS=0; 
        DELS=0;
    end
    dRdT2 = (INE+PINF+PGO+VIS+SUR+ELS+DELS+DPA)/DEN;
    dR_dt=[ r2; dRdT2 ];
end
 
function dpa = pa(t,f,A)
% PRESSURE DIFFERENCE CALCULATION
% LETTY YOU CHANGE THIS TO WHAT YOU NEED TO REPLICATE RESULTS IN PAPERS
%   dt = 1E-6;
%   om = f*2*pi;
%     mn = 3.7;
%     if t < dt - pi/om
%         dpa = 0;
%     elseif t > dt + pi/om
%         dpa = 0;
%     else
%         dpa = -(A*mn*om*sin(om*(dt - t))*...
%             (cos(om*(dt - t))/2 + 1/2)^(mn - 1))/2;
%     end
 
dpa = 0;
 
end
 
function p = pinf(t,f,A)
% FORCING PRESSURE CALCULATION
% LETTY YOU CHANGE THIS TO WHAT YOU NEED TO REPLICATE RESULTS IN PAPERS
%   dt = 1E-6;
%   om = f*2*pi;
%   mn = 3.7;
%     if t < dt - pi/om
%         p = 0;
%     elseif t > dt + pi/om
%         p = 0;
%     else
%         p = -A*(0.5 + 0.5*cos(om*(t - dt))).^mn;
%     end
     
    p = A;
     
end