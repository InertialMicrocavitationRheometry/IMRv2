function [t2,R2,p2,U2,U2dot,T2] = kmsolver1(TFin,Rref,R0,rho8,c8,pAmp,...
    pgo,rpe,enth,gil,neohook,voigt,poly,muval,Gval,Sval,kappa,azero,bzero)
%UNTITLED Summary of this function goes here
% waveform parameters
omega = 0;
pV = 2300;
p8 = 101325;
pA = pAmp*p8;
S = Sval;
TW = 0; % gaussian width (s)
DT = 0; % delay (s)
% viscoelastic properties
mu = muval;%1e-3;
G = Gval;
% yangc = 1;
lambda1 = 0;
lambda2 = 0;
% computational settings
% spectral = 1;
cold = 0; 
vapor = 1;
Nt = 20; Mt = 30; Lt = 3;
Nv = 20; Lv = 3; 
divisions = 1000; 
detail = 5000;
pdisp = 0;
method = '45';
% voigt = 1;
pb0 = pgo*(Rref/R0)^(3*kappa)+pV+2*S/R0;
U0 = sqrt(pb0/rho8);
a0 = azero;
b0 = bzero;

% [t2,R2,p2,U2,~,~,~,T2,~,~,~,~,~,~,~,U2dot] = spectralkm1('detail',detail,'method',method,'yangc',yangc,...
%         'dimout',1,'pdisp',0,'divisions',divisions,... % output
%         'poly',poly,'cold',cold,'vapor',vapor,'nt',Nt,'rpe',rpe,... % secondary effects
%         'tfin',TFin,'pa',pA,'omega',omega,'tw',TW,'DT',DT,... % forcing
%         'pV',pV,'rref',Rref,'r0',R0,'u0',U0,'p0',pgo,'a0',a0,'b0',b0,'p8',p8,'rho8',rho8,'c8',c8,'kappa',kappa,'surf',S,... % physical constants
%         'mu',mu,'g',G,'lambda1',lambda1,'lambda2',lambda2); % viscoelastic properties   

[t2,R2,p2,U2,~,~,~,T2,~,~,~,~,~,~,~,U2dot] = spectralkm2('detail',detail,'method',method,...
        'voigt',voigt,'neohook',neohook,...
        'dimout',1,'pdisp',pdisp,'divisions',divisions,... % output
        'poly',poly,'cold',cold,'vapor',vapor,'nt',Nt,...% secondary effects
        'rpe',rpe,'enth',enth,'gil',gil,... % bubble wall model
        'enth',enth,'tfin',TFin,'pa',pA,'omega',omega,'tw',TW,'DT',DT,... % forcing
        'pV',pV,'rref',Rref,'r0',R0,'u0',U0,'p0',pgo,'a0',a0,'b0',b0,... % physical constants
        'p8',p8,'rho8',rho8,'c8',c8,'kappa',kappa,'surf',S,... % physical constants
        'mu',mu,'g',G,'lambda1',lambda1,'lambda2',lambda2,'plot',0); % viscoelastic properties   
end