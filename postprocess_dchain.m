function [tsol_new, rsol_new, rdotsol_new] = postprocess_dchain(R_eqd,A,w,We,Ca,rho,po,pv,c,tfinal,rmodel,emodel,vmodel,vmaterial)

    
    tspan = [0 tfinal];
    tsegments = 3000;
    tstart = tspan(1);
    tend = tspan(2)/tsegments;
    t = [];
    R = [];  
    
    c2=c^2;                     
    k=1;                        
    S =(rho*c2*R_eqd)/We;       
    eta=(rho*c^2)/Ca;          
    pGo=po-pv+(2*S)/R_eqd;      
    xi = 1/c; 
    r_dt=(pv-po)/(rho*c);
    r_initial = [R_eqd r_dt];
    options = odeset('MaxStep', tfinal/1000);
    
    for i = 1:tsegments
        
        tspansegment = [tstart, tend];
    
        
        [t_add, R_add] = ode23tb(@(t,r) keller_miksis(t,eta,r,k,xi,rho,R_eqd,...
            S,A,w,po,pv,pGo,c,rmodel,emodel,vmodel,vmaterial),tspansegment,r_initial,options);
        
        if(R_add(2,1)*R_add(2,end) < 0)
            
           options = odeset('MaxStep', 1E-9);
           [t_add, R_add] = ode23tb(@(t,r) keller_miksis(t,eta,r,k,xi,rho,R_eqd,...
            S,A,w,po,pv,pGo,c,rmodel,emodel,vmodel,vmaterial),tspansegment,r_initial,options);
        else
            options = odeset('MaxStep', tfinal/1000);   
        end
        
        t = [t t_add'];
        R = [R R_add'];
        r_initial = [R_add(1,end), R_add(2,end)];
        

    end

    tsol_new = t;
    rsol_new = R(1,:);
    rdotsol_new = R(2,:);

end

function [ dR_dt ] = keller_miksis(t,eta,r,k,xi,rho,Req,S,A,w,...
    po,pv,pGo,c,rmodel,emodel,vmodel,vmaterial)
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