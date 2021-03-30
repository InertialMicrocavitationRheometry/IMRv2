function [E0,EPD,LKE,LPE,BIE,VE,TE] = f_energy(Tout,Rout,Pinf,DPinf,Sd,rho,...
    pGo,kappa,po,vmaterial,force)
%F_ENERGY Summary of this function goes here
%   Detailed explanation goes here

%Radius calculations
R = Rout(:,1);
Rdot = Rout(:,2);
R0 = R(1,1);
V  = (4*pi/3)*R.^3;
V0 = (4*pi/3)*R0^3;
PG = pGo*(R0./R).^(3*kappa);

if strcmp(force,'mono')==1
    EPD = 0;
elseif strcmp(force,'sine')==1
    EPD = cumtrapz(Tout,DPinf.*V);
end

LKE = 2*pi*rho*R.^3.*Rdot.^2;
LPE = Pinf.*V;
BIE = PG.*V/(kappa-1)+Sd*(4*pi*R.^2);

E0 = (Pinf(1)+pGo/(kappa-1))*V0+Sd*(4*pi*R0.^2);

[~,mu_inf,mu_o,nc,lambda] = f_carreau(vmaterial,1);
if (strcmp('mu_inf',vmaterial) == 1 || strcmp('mu_0',vmaterial) == 1)
    V = zeros(size(dt));
else
    a = (nc-1)/2;
    V = (mu_inf + (mu_o-mu_inf)*(1+4*lambda^2*(Rdot./R).^2).^a).*Rdot.^2.*R;
end
VE = 16*pi*cumtrapz(Tout,V);

TE  = LKE+LPE+BIE+VE-EPD;
end