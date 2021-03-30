clear; close all; clc;
%INPUT PARAMETERS
Ro_w = 0.5E-6;                      % initial bubble radius = 1E-6  
Sd  = 0.075;                        % surface tension                
po  = 1E5;                          % atmospheric pressure = 101325
pv  = 2300;                         % pressure of vapor          
rho = 1000;                         % density of liquid = 1000      
c   = 1500;                         % speed of sound 
kappa = 1.4;                        % ratio of specific heats
pGo=po-pv+(2*Sd)/Ro_w;
fnatural = (1/(2*pi*Ro_w))*sqrt((3*kappa*(po-pv)+...
    (3*kappa-1)*(2*Sd)/Ro_w)/rho);
tRC = 1/fnatural;                 % natural period
tmag = 20.01;
tfinal=tmag*tRC;                % final time of simulation
deltap = 0.25E6;                   % amplitude of wave
f = 2E6; %1/tRC;                 % driving frequecy
shearmodulus = 0;               % shear modulus of the surounding material
mu_water = 8.9*10E-4;           % viscosity of water
Ca = (rho*c*c)/shearmodulus;    % shear modulus of soft material 
We_w = (rho*c*c*Ro_w)/Sd;       % Weber number, surface tension

%SETTING UP FIGURES
p_figsetup;

% SIMULATION SETTINGS
rmodel = 'RP';
emodel = 'NeoH';
force = 'sine';
filesuffix = '_baseline';
%RUNNING AND PLOTTING
% blood carreau
rstarlim = 1000;
rstarres = 1;
rstarlimtick = 10;
mufilter = 1;
vmaterial = 'lsq_blood';
contourshift = 1;
[Tout,Rout,Rddot,Pinf]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmaterial);
lm = 'k-';
p_figplot;
p_ffield;
% p_figenergy;

%Save figures
fn = strcat('./bubblewallfigs/fRofT',filesuffix);
figure(1)
saveas(gcf,fn,'png')

fn = strcat('./bubblewallfigs/fgammaofT',filesuffix);
figure(2)
saveas(gcf,fn,'png')

fn = strcat('./bubblewallfigs/fmuofT',filesuffix);
figure(3)
saveas(gcf,fn,'png') 

fn = strcat('./bubblewallfigs/ftauofT',filesuffix);
figure(4)
saveas(gcf,fn,'png')

fn = strcat('./bubblewallfigs/ftauofgamma',filesuffix);
figure(5)
saveas(gcf,fn,'png')

fn = strcat('./contourfigs/fcmuofr_muc',filesuffix);
figure(6+contourshift)
saveas(gcf,fn,'png')
fn = strcat('./contourfigs/fcshearofr_muc',filesuffix);
figure(7+contourshift)
saveas(gcf,fn,'png')
fn = strcat('./contourfigs/fctauofr_muc',filesuffix);
figure(8+contourshift)
saveas(gcf,fn,'png')