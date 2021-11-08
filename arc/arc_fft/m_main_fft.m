clear; close all; clc;
routines = strcat(pwd,'/routines');
addpath(routines);
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
tmag = 250.01;
tfinal=tmag*tRC;                % final time of simulation
deltap = 1E3;                     % amplitude of wave
f = 2E6; %1/tRC;                % driving frequecy
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
filesuffix = '_fft';
%RUNNING AND PLOTTING
% blood carreau
rstarlim = 1000;
rstarres = 1;
rstarlimtick = 10;
mufilter = 1;
vmaterial = 'lsq_blood';
contourshift = 1;

[Tout,Rout,~,~]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmaterial);
lm = 'k-';
% p_figplot;
figure(6)
hold on;
f_fourier(Tout,Rout,Ro_w,'-k');

vmaterial = 'mu_0';
[Tout,Rout,~,~]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmaterial);
lm = 'r--';
p_figplot;
figure(6)
f_fourier(Tout,Rout,Ro_w,'r--');

vmaterial = 'mu_inf';
[Tout,Rout,~,~]=f_bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,...
    c,tfinal,force,rmodel,emodel,vmaterial);
figure(6)
f_fourier(Tout,Rout,Ro_w,'-.g');
lm = '-.g';
p_figplot;

figure(6)
set(gca,'XScale','log')
% set(gca,'YScale','log')
xlabel('$f$ [Hz]', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$A$', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
xlim([1E4 1E8])
box on;

%Save figures
fn = strcat('./bubblewallfigs/fRofT',filesuffix);
figure(1)
saveas(gcf,fn,'png')

fn = strcat('./bubblewallfigs/fft',filesuffix);
figure(6)
saveas(gcf,fn,'png')

rmpath(routines);