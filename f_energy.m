function [E0,LKE,LPE,BIE,VE,TE] = f_energy(Tout,Rout,Pinf,Sd,rho,pGo,kappa,po,vmaterial,force,f,deltap)
%F_ENERGY Summary of this function goes here
%   Detailed explanation goes here

%Radius calculations
R = Rout(:,1);
Rdot = Rout(:,2);
R0 = R(1,1);
dt = [diff(Tout);0];
V  = (4*pi/3)*R.^3;
V0 = (4*pi/3)*R0^3;
PG = pGo*(R0./R).^(3*kappa);

if strcmp(force,'mono')==1
    EPD = 0;
elseif strcmp(force,'sine')==1
    Dp = 2*pi*f*deltap*cos(2*pi*f*Tout);
    EPD = cumtrapz(Dp.*V.*dt);
end

LKE = 2*pi*rho*R.^3.*Rdot.^2;
LPE = Pinf.*V;
% LPE = 4*pi*Pinf.*Rdot.*R.^2
BIE = PG.*V/(kappa-1)+Sd*(4*pi*R.^2);

E0 = Pinf(1)*V0+pGo*V0/(kappa-1)+Sd*(4*pi*R0.^2)+EPD;

ve = Rdot.^2.*R;
[~,mu_inf,mu_o,nc,lambda] = f_carreau(vmaterial,1);
if (strcmp('mu_inf',vmaterial) == 1 || strcmp('mu_0',vmaterial) == 1)
    S = zeros(size(dt));
else
    a = (nc-1)/2;
    V = (mu_inf + (mu_o-mu_inf)*(1+4*lambda^2*(Rdot./R).^2).^a).*Rdot.^2.*R;
    S = V.*dt;
end
VE = 16*pi*cumtrapz(S);

TE  = LKE+LPE+BIE+VE;
end