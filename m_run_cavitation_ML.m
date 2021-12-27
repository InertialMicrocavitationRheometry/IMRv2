% This script can be used as a template to run RP_Cav
clear; close all; clc;
% Adding the routines path
routines = strcat(pwd,'/src');
post = strcat(pwd,'/post_processing');
addpath(genpath(routines));
addpath(post);

% All quantities below are in SI units (s,m)
% time to run simulation 
tend = 0.75E-6;
% Initital Radii 
R0 = 0.5E-6;
% Ammount of nodes inside the bubble (~100 is a good to start)
NT = 5; 
% Ammount of nodes outside the bubble (~100 is a good to start)
NTM = 5; 
% Type of external forcing. Refer to RP_Cav
Pext_type = 'ML'; 
% [ Pressure ; Freq ] 
Pext_Amp_Freq =[5E5 0]; 
% 1 : display simulation time, 0 : do not display
disptime = 1; 
% Thermal effects inside bubble, 1: yes, 0: no
Tgrad = 0; 
% Thermal effects outside bubble, 1: yes, 0: no
Tmgrad = 0;
% Vapor diffusion effects, 1: yes, 0: no
Cgrad = 0;  
% Output variables in dimensional form, 1: yes, 0: no
Dim = 0; 
% Activates the effect of compressibility 
% 1: Keller-Miksis w/ pressure, 0: Rayleigh-Plesset
comp = 1; 
% material to calculate viscosity
vmaterial = 'blood_combined';
% non-Newtonian model for viscosity
vmodel = 'carreau';

% RUNNING BASELINE CODE
[ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);

% POST PROCESSING RESULTS
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27);
lm = 'k';
r_fig_ML;
r_ffield_ML;
rmpath(routines);
rmpath(post);