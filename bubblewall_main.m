clear; close all; clc;
%INPUT PARAMETERS
Ro_w=1E-4;                              % initial bubble radius = 1E-6  
Sd=0.072;                               % surface tension                
po=101325;                              % atmospheric pressure = 101325
pv=2300;                                % pressure of vapor          
rho=1060;                               % density of liquid = 1000      
mud=1E-3;                               % dynamic viscosity
c=1500;                                 % speed of sound 
A=0.0001*po;%15E6;                       % amplitude of wave
w=345E3;                                % frequency of wave = 345E32
tRC = 0.915*sqrt(rho*Ro_w^2/(po+A-pv)); % Taken from Brennen equation 2.40
tfinal=4*tRC;                           % final time of simulation
shearmodulus = 0;                       % shear modulus of the surounding material
mu_water = 8.9*10E-4;                   % viscosity of water
Ca = (rho*c*c)/shearmodulus;            % shear modulus of soft material 
We_w = (rho*c*c*Ro_w)/Sd;
gamma0 = A/(rho*c*Ro_w);

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

figure(5)
hold on
xlabel('$\dot{\gamma}$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_v$ [Pa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')


%RUNNING AND PLOTTING
%Water
% [ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','water'); 
% plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'b','LineWidth',2); 
%blood inf

mu_o = 0.112; 

vmaterial = 'lsq_mu_inf';
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau',vmaterial);
figure(1)
tspan = (ODE.x/tRC);
Roft = (ODE.y(1,:)/Ro_w)-1;
plot(tspan,Roft, 'g','LineWidth',2); 
figure(2)
gamma = (ODE.y(2,:))./(ODE.y(1,:));
gamman = gamma/gamma0;
plot(tspan,gamman,'g','LineWidth',2); 
figure(3)
mutrace = carreau(vmaterial,gamma)/mu_o*ones(size(gamma));
plot(tspan,mutrace,'g','LineWidth',2); 
figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace,'g','LineWidth',2); 
figure(5)
tautrace = abs(mutrace.*gamma).*(gamma>0);
plot(gamma,tautrace,'g','LineWidth',2); 

%Blood
vmaterial = 'lsq_blood';
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau',vmaterial); %blood - Carreau 
figure(1)
tspan = (ODE.x/tRC);
Roft = (ODE.y(1,:)/Ro_w)-1;
plot(tspan,Roft,'r--','LineWidth',2); 
figure(2)
gamma = (ODE.y(2,:))./(ODE.y(1,:));
gamman = gamma/gamma0;
plot(tspan,gamman,'r--','LineWidth',2);
figure(3)
mutrace = carreau(vmaterial,gamma)/mu_o;
plot(tspan,mutrace, 'r--','LineWidth',2);
figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace, 'r--','LineWidth',2);
figure(5)
tautrace = abs(mutrace.*gamma).*(gamma>0);
plot(gamma,tautrace, 'r--','LineWidth',2);

%blood 0
vmaterial = 'lsq_mu_knot';
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau',vmaterial);
figure(1)
tspan = (ODE.x/tRC);
Roft = (ODE.y(1,:)/Ro_w)-1;
plot(tspan,Roft,'k-.','LineWidth',2); 
figure(2)
gamma = (ODE.y(2,:))./(ODE.y(1,:));
gamman = gamma/gamma0;
plot(tspan,gamman,'k-.','LineWidth',2); 
figure(3)
mutrace = carreau(vmaterial,gamma)/mu_o*ones(size(gamma));
plot(tspan,mutrace,'k-.','LineWidth',2); 
figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace,'k-.','LineWidth',2); 
figure(5)
tautrace = abs(mutrace.*gamma).*(gamma>0);
plot(gamma,tautrace,'k-.','LineWidth',2); 
xlim([0 max(gamma)])

%Save figures
figure(1)
print('fRofT', '-dpng'); 
figure(2)
print('fgammaofT', '-dpng'); 
figure(3)
print('fmuofT', '-dpng'); 
figure(4)
print('ftauofT', '-dpng'); 
figure(5)
print('ftauofgamma', '-dpng'); 