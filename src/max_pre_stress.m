function [Smaxpred] = max_pre_stress(Ro, al_nd, pwv_nd, We, Re, De, Ca, alpha)
    trc = sqrt(6*pi)*gamma(11/6)/(5*gamma(4/3));
    fbarst = pi./(sqrt(6).*We.*trc);
    gam = 1.4;
    B = 2.1844;
    fbarbc = -(1-pwv_nd+1./(We.*Ro)).*B.*Ro.^(3*gam)+pwv_nd;
    Mc = 1./al_nd;
    fbarc = -2.*Mc./(Mc+sqrt(Mc.^2+4.*trc^2));
    C = 0.46379+0.56391./Re+5.74916./Re.^2;
    fbarv = -4.*C.^2./(2.*C.^2+sqrt(4.*C.^4+C.^2.*Re.^2.*trc^2));
    fbarmax = fbarv+De./trc.*((fbarv).*exp(-trc./De) -fbarv);
    fbare = 1/(60*Ca*gamma(5/6))*gamma(1/3)*((40*sqrt(pi)*Ro*(1-3*alpha))+ ...
        (120*(-1+2*Ro^3)*alpha*gamma(7/6))/(Ro*gamma(2/3))+((-50+177*alpha))*gamma(5/6)/gamma(4/3));
    fbarsls = fbarmax-fbare;
    fsum =fbarbc+fbarst+fbarc+fbarsls;
    tg = (5*sqrt(pi)*gamma(5/6)-6*Ro^(5/2)*gamma(4/3)*hypergeom([1/2,5/6],11/6,Ro^3))./(5*sqrt(6-6.*fsum)*gamma(4/3));
    Smaxpred = fbarv.*(1-exp(-tg./De));
end