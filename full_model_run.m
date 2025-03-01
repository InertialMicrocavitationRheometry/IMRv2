clc;
clear;
close;

addpath('src');

% equation options
R0 = 2.447495043190468e-04;
Req = 3.008409399929516e-05;
% Req = 0.122917895573930*R0;
ratio = Req/R0;
% Req = R0;
tfin = R0/3;
% pzero = 0.001887748193919;
% cstar = 1.520710024149080e+02;
kappa = 1.4;
Lheat = 2.378193575129533e+04;
T8 = 298.15;
rho8 = 1064;
tvector = linspace(0,tfin,1000);
radial = 2;
vapor = 1;
bubtherm = 1;
medtherm = 1;
masstrans = 1;
stress = 1;
varin = {'progdisplay',0,...
    'radial',radial,...
         'bubtherm',bubtherm,...
         'tvector',tvector,...
         'vapor',vapor,...
         'medtherm',medtherm,...
         'method',23,...
         'stress',stress,...
         'mu',0.027606000000000,...
         'g',3.125000000000000e+02,...
         'lambda1',0,...
         'r0',R0,...
         'req',Req,...
         'kappa',kappa,...
         're',92.055015921905493,...
         'masstrans',masstrans,...
         't8',T8,...
         'rho8',rho8};

[tf,Rf,~] = m_imr_full_model(varin{:},'Nt',50,'Mt',50);
% [td,Rd,~] = m_imr_finitediff(varin{:},'Nt',50,'Mt',50);

figure(1)
hold on;
plot(tf,Rf,'^')
% plot(td,Rd,'s')
ylim([0 1.2])