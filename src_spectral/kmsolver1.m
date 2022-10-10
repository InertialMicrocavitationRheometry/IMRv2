% function [t2,R2,p2,U2,U2dot,T2] = kmsolver1(TFin,Rref,R0,rho8,c8,pAmp,...
%     pgo,rpe,enth,gil,neohook,voigt,poly,muval,Gval,Sval,kappa,azero,bzero)
%UNTITLED Summary of this function goes here
% waveform parameters
clear; close all; clc;

TFin = 1E-6;
Rref = 1E-6;
R0 = 1E-6;
options = optimset('TolX',1E-10,'MaxIter',50);
rho8 = 1000;
kappa = 1.4;
c8 = 1500;
pgo = 2300;
%stress model
voigt = 0;
neohook = 1;
%bubble wall model
rpe = 0;
enth = 1;
gil = 0;
poly = 1;
Sval = 0.07;
muval = 0.115;
% muval = 1E-3;
Gval = 0;
alpha0 = 1; 
epsiloncrit = 0.29;

omega = 0;
pV = 2300;
p8 = 101325;
pA = 1*p8;
p0 = p8 + 2*Sval/Rref - pV;
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
U0 = -sqrt((pA+p8)/rho8);
a0 = [];
b0 = [];

% [t2,R2,p2,U2,~,~,~,T2,~,~,~,~,~,~,~,U2dot] = spectralkm1('detail',detail,'method',method,'yangc',yangc,...
%         'dimout',1,'pdisp',0,'divisions',divisions,... % output
%         'poly',poly,'cold',cold,'vapor',vapor,'nt',Nt,'rpe',rpe,... % secondary effects
%         'tfin',TFin,'pa',pA,'omega',omega,'tw',TW,'DT',DT,... % forcing
%         'pV',pV,'rref',Rref,'r0',R0,'u0',U0,'p0',pgo,'a0',a0,'b0',b0,'p8',p8,'rho8',rho8,'c8',c8,'kappa',kappa,'surf',S,... % physical constants
%         'mu',mu,'g',G,'lambda1',lambda1,'lambda2',lambda2); % viscoelastic properties   

[t2,R2,p2,U2,~,~,~,T2,~,~,~,~,~,~,~] = f_imrv2('detail',detail,'method',method,...
        'voigt',voigt,'neohook',neohook,...
        'dimout',1,'pdisp',pdisp,'divisions',divisions,... % output
        'poly',poly,'cold',cold,'vapor',vapor,'nt',Nt,...% secondary effects
        'rpe',rpe,'enth',enth,'gil',gil,... % bubble wall model
        'enth',enth,'tfin',TFin,'pa',pA,'omega',omega,'tw',TW,'DT',DT,... % forcing
        'pV',pV,'rref',Rref,'r0',R0,'u0',U0,'p0',p0,'a0',a0,'b0',b0,... % physical constants
        'p8',p8,'rho8',rho8,'c8',c8,'kappa',kappa,'surf',S,... % physical constants
        'mu',mu,'g',G,'lambda1',lambda1,'lambda2',lambda2,'plot',0); % viscoelastic properties   
    
plot(t2/TFin,R2/R0)
% end