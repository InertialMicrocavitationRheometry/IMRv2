% This script can be used as a template to run RP_Cav
clear all; close all; clc;
% Adding the routines path
routines = strcat(pwd,'/src');
addpath(routines);

% All quantities below are in SI units (s,m)
% time to run simulation 
tspan = 5e-3;                   
% Initital Radii 
R0 = 500e-6;                    
% Ammount of nodes inside the bubble (~100 is a good to start)
NT = 100; 
% Ammount of nodes outside the bubble (~100 is a good to start)
NTM = 120; 
% Type of external forcing. Refer to RP_Cav
Pext_type = 'GS'; 
% [ Pressure ; Freq ] 
Pext_Amp_Freq =[101325 1E6]; 
% 1 : display simulation time, 0 : do not display
disptime = 1; 
% Thermal effects inside bubble, 1: yes, 0: no
Tgrad = 1; 
% Thermal effects outside bubble, 1: yes, 0: no
Tmgrad = 1;
% Vapor diffusion effects, 1: yes, 0: no
Cgrad = 1;  
% Output variables in dimensional form, 1: yes, 0: no
Dim = 1; 
% Activates the effect of compressibility, 1: Keller-Miksis w/ pressure, 0:
% Rayleigh-Plesset
comp = 1; 
% material to calculate viscosity
vmaterial = 'water';
% non-Newtonian model for viscosity
vmodel = 'newtonian';

[ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tspan,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);

figure(1)
hold on;
plot(t,R,'k')
figure(2)
plot(t,P,'k')
rmpath(routines);
