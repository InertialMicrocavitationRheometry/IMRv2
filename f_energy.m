function [E0,LKE,LPE,BIE,VE,TE] = f_energy(Tout,Rout,Pinf,Sd,rho,pGo,kappa,po,vmaterial)
%F_ENERGY Summary of this function goes here
%   Detailed explanation goes here

%Radius calculations
R = Rout(:,1);
Rdot = Rout(:,2);
R0 = R(1,1);

V  = (4*pi/3)*R.^3;
V0 = (4*pi/3)*R0^3;
PG = pGo*(R0./R).^(3*kappa);

E0 = Pinf(1)*V0+pGo*V0/(kappa-1)+Sd*(4*pi*R0.^2);

LKE = 2*pi*rho*R.^3.*Rdot.^2;
LPE = Pinf.*V;
BIE = PG.*V/(kappa-1)+Sd*(4*pi*R.^2);
ve = Rdot.^2.*R;
dt = [0;diff(Tout)];

[mu,mu_inf,mu_o,nc,lambda] = f_carreau(vmaterial,1);

if (strcmp('mu_inf',vmaterial) == 1 || strcmp('mu_0',vmaterial) == 1)
    S = zeros(size(dt));
else
    S = R.^3.*ve.*f_visG(nc,lambda,Rdot,R).*dt;
    mu = mu_inf;
end

VE = 16*pi*mu.*cumtrapz(ve.*dt) + 48*pi*(mu_o-mu_inf)*cumtrapz(S);
TE  = LKE+LPE+BIE+VE;
end