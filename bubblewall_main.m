clear; close all; clc;
%INPUT PARAMETERS
Ro_w=10E-6;                              % initial bubble radius = 1E-6  
Sd=0.072;                               % surface tension                
po=101325;                              % atmospheric pressure = 101325
pv=2300;                                % pressure of vapor          
rho=1060;                               % density of liquid = 1000      
mud=1E-3;                               % dynamic viscosity
c=1500;                                 % speed of sound 
kappa = 1.4;                            % ratio of specific heats
omegan = (1/(2*pi*Ro_w))*sqrt((3*kappa*(po-pv)+(3*kappa-1)*(2*Sd)/Ro_w)/rho);
tRC = 1/omegan;                         % natural period
tfinal=1*tRC;                          % final time of simulation
deltap = 1E0;                           % amplitude of wave
f = (2)*1/tRC;                         % driving frequecy
shearmodulus = 0;                       % shear modulus of the surounding material
mu_water = 8.9*10E-4;                   % viscosity of water
Ca = (rho*c*c)/shearmodulus;            % shear modulus of soft material 
We_w = (rho*c*c*Ro_w)/Sd;
mu_o = 0.056; 


% SETTING UP THE FIGURES
figure(1)  
hold on
xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{R^*}$-1', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xticks([0:0.25:1])
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;

figure(2)
hold on
xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\dot{\gamma}^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks([0:0.25:1])
% ylim([-1.5 1.5]*1E-5)

figure(3)
hold on
xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\mu$*', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
ylim([0 1.1])
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xticks([0:0.25:1])
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;

figure(4)
hold on
xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_v$ [Pa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xticks([0:0.25:1])
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;

figure(5)
hold on
xlabel('$\dot{\gamma}$ [1/s]', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_v$ [Pa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;

%RUNNING AND PLOTTING
%Water
% [ODE]=bubblewall_solver(Ro_w,deltap,kappa,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','water'); 
% plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'b','LineWidth',2); 
%blood inf

rmodel = 'KM';
emodel = 'NeoH';
vmodel = 'Carreau';
vmaterial = 'mu_inf';
force = 'sine';
[ODE]=bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,c,tfinal,...
    force,rmodel,emodel,vmodel,vmaterial);
figure(1)
tspan = ODE.x*omegan;
Roft = (ODE.y(1,:)/Ro_w)-1;
plot(tspan,Roft, 'g-','LineWidth',2); 
figure(2)
gamma = (ODE.y(2,:))./(ODE.y(1,:));
gamman = gamma*tRC;
plot(tspan,gamman,'g-','LineWidth',2); 
figure(3)
mutrace = carreau(vmaterial,gamma)/mu_o*ones(size(gamma));
plot(tspan,mutrace,'g-','LineWidth',2); 
figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace,'g-','LineWidth',2); 
figure(5)
tautrace = mutrace.*gamma;
taug = abs([gamma;tautrace]);
taug = sortrows(taug');
gt = taug(:,1);
taut = taug(:,2);
plot(gt,taut,'g-','LineWidth',2); 

%Blood
vmaterial = 'lsq_blood';
[ODE]=bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,c,tfinal,...
    force,rmodel,emodel,vmodel,vmaterial);
figure(1)
tspan = ODE.x*omegan;
Roft = (ODE.y(1,:)/Ro_w)-1;
plot(tspan,Roft,'r--','LineWidth',2); 
figure(2)
gamma = (ODE.y(2,:))./(ODE.y(1,:));
gamman = gamma*tRC;
plot(tspan,gamman,'r--','LineWidth',2);
figure(3)
mutrace = carreau(vmaterial,gamma)/mu_o;
plot(tspan,mutrace, 'r--','LineWidth',2);
figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace, 'r--','LineWidth',2);
figure(5)
tautrace = mutrace.*gamma;
taug = abs([gamma;tautrace]);
taug = sortrows(taug');
gt = taug(:,1);
taut = taug(:,2);
plot(gt,taut, 'r--','LineWidth',2);

%blood 0
vmaterial = 'mu_knot';
[ODE]=bubblewall_solver(Ro_w,deltap,kappa,f,We_w,Ca,rho,po,pv,c,tfinal,...
    force,rmodel,emodel,vmodel,vmaterial);
figure(1)
tspan = ODE.x*omegan;
Roft = (ODE.y(1,:)/Ro_w)-1;
plot(tspan,Roft,'k-.','LineWidth',2); 
figure(2)
gamma = (ODE.y(2,:))./(ODE.y(1,:));
gamman = gamma*tRC;
plot(tspan,gamman,'k-.','LineWidth',2); 
figure(3)
mutrace = carreau(vmaterial,gamma)/mu_o*ones(size(gamma));
plot(tspan,mutrace,'k-.','LineWidth',2); 
figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace,'k-.','LineWidth',2); 
figure(5)
tautrace = mutrace.*gamma;
taug = abs([gamma;tautrace]);
taug = sortrows(taug');
gt = taug(:,1);
taut = taug(:,2);
plot(gt,taut,'k-.','LineWidth',2); 
% xlim([0 max(gamma)])

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
