clear; close all; clc;
%INPUT PARAMETERS
Ro_w = 0.5E-6;                      % initial bubble radius = 1E-6  
Sd  = 0.072;                       % surface tension                
po  = 1E5;                      % atmospheric pressure = 101325
pv  = 2300;                        % pressure of vapor          
rho = 1000;                       % density of liquid = 1000      
c   = 1500;                         % speed of sound 
kappa = 1.4;                    % ratio of specific heats
fnatural = (1/(2*pi*Ro_w))*sqrt((3*kappa*(po-pv)+...
    (3*kappa-1)*(2*Sd)/Ro_w)/rho);
tRC = 1/fnatural;                 % natural period
tmag = 10.0001;
tfinal=tmag*tRC;                % final time of simulation
deltap = 0.2E6;                   % amplitude of wave
f = 6.4E6;%1/tRC;                 % driving frequecy
shearmodulus = 0;               % shear modulus of the surounding material
mu_water = 8.9*10E-4;           % viscosity of water
Ca = (rho*c*c)/shearmodulus;    % shear modulus of soft material 
We_w = (rho*c*c*Ro_w)/Sd;       % Weber number, surface tension
mu_o = 0.056;                   % blood no shear viscosity
stress = po;

%SETTING UP FIGURES
figsetup;

% SIMULATION SETTINGS
rmodel = 'RP';
emodel = 'NeoH';
vmodel = 'Carreau';
force = 'sine';

%RUNNING AND PLOTTING
vmaterial = 'mu_inf';
[Tout,Rout,Rddot]=bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmodel,vmaterial);
lm = 'g-';
figplot;

%Blood
vmaterial = 'lsq_blood';
[Tout,Rout,Rddot]=bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmodel,vmaterial);
lm = 'r--';
figplot;
rstarlim = 500;
ffield;

%blood 0
vmaterial = 'mu_0';
[Tout,Rout,Rddot]=bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmodel,vmaterial);
lm = 'k-.';
figplot;

%Save figures
figure(1)
saveas(gcf,'fRofT','png')
figure(2)
saveas(gcf,'fgammaofT','png')
figure(3)
saveas(gcf,'fmuofT','png') 
figure(4)
saveas(gcf,'ftauofT','png')
figure(5)
saveas(gcf,'ftauofgamma','png')
figure(6)
saveas(gcf,'fmuofr','png')
figure(7)
saveas(gcf,'fshearofr','png')
figure(8)
saveas(gcf,'ftauofr','png')