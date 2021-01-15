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
mu_inf = 0.00345;
stress = po;

%SETTING UP FIGURES
m_figsetup;

% SIMULATION SETTINGS
rmodel = 'RP';
emodel = 'NeoH';
vmodel = 'Carreau';
force = 'sine';

%RUNNING AND PLOTTING

%blood inf
rstarlim = 5;
rstarres = 80;
rstarlimtick = rstarlim/5;
vmaterial = 'mu_inf';
mufilter = 0;
contourshift = 0;
[Tout,Rout,Rddot]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmodel,vmaterial);
lm = 'g-';
m_figplot;
m_ffield;

% blood carreau
rstarlim = 700;
rstarres = 5;
rstarlimtick = 100;
mufilter = 1;
vmaterial = 'lsq_blood';
contourshift = 3;
[Tout,Rout,Rddot]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmodel,vmaterial);
lm = 'r--';
m_figplot;
m_ffield;

% blood o
rstarlim = 5;
rstarres = 80;
rstarlimtick = rstarlim/5;
vmaterial = 'mu_0';
mufilter = 0;
contourshift = 6;
[Tout,Rout,Rddot]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmodel,vmaterial);
lm = 'k-.';
m_figplot;
m_ffield;

%Save figures
figure(1)
saveas(gcf,'./bubblewallfigs/fRofT','png')
figure(2)
saveas(gcf,'./bubblewallfigs/fgammaofT','png')
figure(3)
saveas(gcf,'./bubblewallfigs/fmuofT','png') 
figure(4)
saveas(gcf,'./bubblewallfigs/ftauofT','png')
figure(5)
saveas(gcf,'./bubblewallfigs/ftauofgamma','png')

figure(7)
saveas(gcf,'./contourfigs/fcshearofr_mui','png')
figure(8)
saveas(gcf,'./contourfigs/fctauofr_mui','png')
figure(9)
saveas(gcf,'./contourfigs/fcmuofr_muc','png')
figure(10)
saveas(gcf,'./contourfigs/fcshearofr_muc','png')
figure(11)
saveas(gcf,'./contourfigs/fctauofr_muc','png')
figure(13)
saveas(gcf,'./contourfigs/fcshearofr_muo','png')
figure(14)
saveas(gcf,'./contourfigs/fctauofr_muo','png')