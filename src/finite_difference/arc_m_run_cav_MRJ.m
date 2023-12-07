% This script can be used as a template to run RP_Cav
clear; close all; clc;
% Adding the routines path
% routines = strcat(pwd,'/src');
% post = strcat(pwd,'/post_processing');
% addpath(genpath(routines));
% addpath(post);

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
Pext_Amp_Freq =[410 0]; 
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

% RUNNING BASELINE CODE
[ t , R ,U ,P, T, C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);

plot(t,R)
%%
% POST PROCESSING RESULTS
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27); Uc = Pmt(28);
lm = 'k';
r_fig_baseline;
r_ffield_baseline;

lm = '-.r';
vmaterial = 'blood_zero';
vmodel = 'newtonian';
[ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27); Uc = Pmt(28);
r_fig_baseline;

lm = '--b';
vmaterial = 'blood_infinity';
vmodel = 'newtonian';
[ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27); Uc = Pmt(28);
r_fig_baseline;

rmpath(routines);
rmpath(post);