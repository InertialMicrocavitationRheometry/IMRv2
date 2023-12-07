% This script can be used as a template to run RP_Cav
clear; close all; clc;
% Adding the routines path
src= strcat(pwd,'/src/finite_difference');
addpath(genpath(src));

% All quantities below are in SI units (s,m)
% time to run simulation 
tend = 1.738E-4; % WARNING: Does not matter, gets changed later
% Initital Radii 
R0 = 3.5E-4;
% Ammount of nodes inside the bubble (~100 is a good to start)
NT = 40; 
% Ammount of nodes outside the bubble (~100 is a good to start)
NTM = 5; 
% Type of external forcing. Refer to RP_Cav
Pext_type = 'IC'; 
% [ Pressure ; Freq ] 
Pext_Amp_Freq =[409 0]; 
% 1 : display simulation time, 0 : do not display
disptime = 0; 
% Thermal effects inside bubble, 1: yes, 0: no
Tgrad = 1; 
% Thermal effects outside bubble, 1: yes, 0: no
Tmgrad = 0;
% Vapor diffusion effects, 1: yes, 0: no
Cgrad = 1;  
% Output variables in dimensional form, 1: yes, 0: no
Dim = 0; 
% Activates the effect of compressibility 
% 1: Keller-Miksis w/ pressure, 0: Rayleigh-Plesset
comp = 1; 
% material to calculate viscosity
vmaterial = 'water';
% non-Newtonian model for viscosity
vmodel = 'newtonian';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% RUNNING BASELINE CODE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,R,~,~,~,~,~,~,~,~] = m_cavitation...
    (tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
    Dim,comp,vmaterial,vmodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [t,R];
%filename = strcat('../../Simulation_Data_',num2str(mu_best),'_',num2str(G_best),'.csv');
%writematrix(a, filename);

figure(1);
plot(t,R,'k')
xlabel('t')
ylabel('R')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Clean up %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath(src);