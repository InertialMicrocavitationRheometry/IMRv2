clear; close; clc;
%INPUT PARAMETERS
Ro_w=10E-6;                             %initial bubble radius = 1E-6  CHECKED
Sd=0.072;                               %surface tension                CHECKED
po=101325;                              %atmospheric pressure = 101325
pv=2300;                                %pressure of vapor          CHECKED
rho=1060;                               %density of liquid = 1000      CHECKED
mud=1E-3;                               %dynamic viscosity
c=1500;                                 %speed of sound CHECKED
A=(1+1E-1)*po;%15E6;                          %amplitude of wave
w=345E3;                                %frequency of wave = 345E32
tRC = 0.915*sqrt(rho*Ro_w^2/(A-po-pv)); % Taken from Brennen equation 2.40
tfinal=200*tRC;                          %final time of simulation
shearmodulus = 0;
mu_water = 8.9*10E-4;                   %viscosity of water
Ca = (rho*c*c)/shearmodulus;            %shear modulus of soft material = ((.1E6)/3 CHECKED
We_w = (rho*c*c*Ro_w)/Sd;


%%
close;
figure(12)  %FIGURE 1 with LEAST SQUARES FIT VALUES
hold on 
%Carreau model for blood found with least squares
mu_o = 0.1117; 
mu_inf = 0.00262; 
lambda = 39.7475;
gammadot = linspace(0,4,1E3);
tau_wo = mu_o.*gammadot;
tau_winf = mu_inf.*gammadot + (mu_o-mu_inf).*(lambda.^(-1/2)).*gammadot.^(1/2);
tau_blood = carreau('lsq_blood',gammadot).*gammadot;
plot(gammadot,tau_blood*1e3,'k','LineWidth',3); 
plot(gammadot,tau_wo*1e3,'b-.','LineWidth',3); 
plot(gammadot,tau_winf*1e3,'r--','LineWidth',3); 
xlabel('$\dot{\gamma} = \dot{R}/R$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\tau_{rr}$ [mPa]','Interpreter','Latex', 'FontSize',20); 
%leg = legend('$\mu_{blood}$','$\mu_{b,0}$','$\mu_{b,\infty}$'); 
%set(leg,'Interpreter','latex','Location','northwest');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f12new2', '-dpng'); 
%%
close;
figure(13)  %FIGURE 2 with LEAST SQUARES FIT VALUES
hold on
tRC = 0.915*sqrt(rho*Ro_w^2/(A-po-pv));   % Taken from Brennen equation 2.40  
%Water
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','water'); 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'b','LineWidth',2); 

%blood inf
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_mu_inf');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'g:','LineWidth',2); 

%Blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood'); %blood - Carreau 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'r--','LineWidth',2); 
%blood 0
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_mu_knot');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'k-.','LineWidth',2); 

xticks([0:5:40]);
yticks([0:.2:1]); 
axis([0 40 0 1]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{R/R$_{0}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
%print('f2New', '-dpng'); 
print('f13New2', '-dpng'); 


%%
close; clc;
figure(14) %FIGURE 4 with LEAST SQUARES FIT VALUES 
hold on 
mu_blood = 0.112;
mu_water = 8.9*10E-4; 
mu_o = 0.112; 
%Carreau model for blood
tfinal = 40*tRC;
% [ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood');


[tsol_new, rsol_new, rdotsol_new, RR] = postprocess_dchain(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood');

plot(tsol_new/tRC, rdotsol_new/Ro_w);
% print('rvstime', '-dpng'); 

figure (100)
plot((tsol_new/tRC), carreau('lsq_blood',rdotsol_new./rsol_new)/mu_blood, '.r','LineWidth',2); 
% line(mu_water/mu_o, 'b','LineWidth',2);
% line(mu_o/mu_o, 'k-.','LineWidth',2)
% line(mu_inf/mu_o, 'g-.','LineWidth',2)
%xticks([0:20:120]);
% yticks([0:.2:1.2]); 
%axis([0 125 0 1.2]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\mu$/$\mu_{0}$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('Carreau'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f14new2', '-dpng'); 

%%
close;
figure(15) %FIGURE 5 with LEAST SQUARES FIT VALUES   
hold on 
%mu_blood = 0.056; 
%Carreau model for blood

[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood'); 
gammadot = rdot_new./r_new;
tau = carreau('lsq_blood',gammadot).*gammadot/1E6;
plot(t_new/tRC,tau,'k','LineWidth',2); 


[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_mu_knot'); 
gammadot = rdot_new./r_new;
tau = carreau('lsq_mu_knot',gammadot).*gammadot/1E6;
plot(t_new/tRC,tau,'--r','LineWidth',2); 

[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_mu_inf'); 
gammadot = rdot_new./r_new;
tau = carreau('lsq_mu_inf',gammadot).*gammadot/1E6;
plot(t_new/tRC,tau,'.g','LineWidth',2); 

%{
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_inf'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('mu_inf',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'b-.','LineWidth',2); 
%}

xticks([0:20:120]);
% yticks([-60]); 
%axis([0 120 -0.5 0.5]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}$ [MPa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{blood}$','$\mu_{blood,0}$', '$\mu_{blood,\infty}$'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f15New2', '-dpng'); 

%%   

%RAYLEIGH PLESSET 


%%

close;
figure(16) %RAYLEIGH PLESSET w/ FIGURE 2
hold on
tRC = 0.915*sqrt(rho*Ro_w^2/(A-po-pv));   % Taken from Brennen equation 2.40  
tfinal=200*tRC;
%Water
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','water'); 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'b','LineWidth',2); 

%blood inf
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_mu_inf');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'g:','LineWidth',2); 

%Blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_blood'); %blood - Carreau 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'r--','LineWidth',2); 
%blood 0
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_mu_knot');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'k-.','LineWidth',2); 

xticks([0:10:70]);
yticks([0:.2:1]); 
axis([0 70 0 1]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{R/R$_{0}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
%print('f2New', '-dpng'); 
print('f16New2', '-dpng'); 


%%
close;
figure(17) %RAYLEIGH PLESSET w/ FIGURE 4 
hold on 
mu_blood = 0.1120;
mu_water = 8.9*10E-4; 
mu_o = 0.112; 
%Carreau model for blood  
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_blood'); 
plot((ODE.x(2:end)/tRC), carreau('lsq_blood',ODE.y(2,2:end)./ODE.y(1,2:end))/mu_blood, 'r','LineWidth',2); 
yline(mu_water/mu_o, 'b','LineWidth',2);
yline(mu_o/mu_o, 'k-.','LineWidth',2)
yline(mu_inf/mu_o, 'g-.','LineWidth',2)
xticks([0:50:200]);
yticks([0:.2:1.2]); 
axis([0 200 0 1.2]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\mu$/$\mu_{0}$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('Carreau'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f17new2', '-dpng'); 

%%
close;
figure(18)  %RAYLEIGH PLESSET w/ FIGURE 5  
hold on 
%mu_blood = 0.056; 
%Carreau model for blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_blood'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('lsq_blood',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'k','LineWidth',2); 


[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_mu_knot'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('lsq_mu_knot',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'--r','LineWidth',2); 


[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'RP','NeoH','Carreau','lsq_mu_inf'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('lsq_mu_inf',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'g-.','LineWidth',2); 


xticks([0:1:5]);
% yticks([-60]); 
axis([0 4 -0.5 0.5]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}$ [MPa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{blood}$','$\mu_{blood,0}$', '$\mu_{blood,\infty}$'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f18New2', '-dpng'); 

%{
%%
close;
figure(19) %FIGURE 15 with longer time
hold on 
%mu_blood = 0.056; 
%Carreau model for blood 
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('lsq_blood',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'k','LineWidth',2); 


[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_mu_knot'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('lsq_mu_knot',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'--r','LineWidth',2); 

[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_mu_inf'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('lsq_mu_inf',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'--g','LineWidth',2); 

%{
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_inf'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('mu_inf',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'b-.','LineWidth',2); 
%}

xticks([0:5:15]);
% yticks([-60]); 
axis([0 15 -0.5 0.5]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}$ [MPa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{blood}$','$\mu_{blood,0}$', '$\mu_{blood,\infty}$'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f19New2', '-dpng'); 

%%
close;
figure(20) %Gamma dot vs time
hold on 


[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
plot((ODE.x/tRC),gammadot, 'b','LineWidth',2); 



%xticks([0:5:15]);
% yticks([-60]); 
%axis([0 15 -0.5 0.5]); 
xlabel('time', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('gammadot', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{blood}$','$\mu_{blood,0}$', '$\mu_{blood,\infty}$'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f20New2', '-dpng'); 

%%
close;
figure(21) %Gamma dot vs time
hold on 

tfinal = 5*tRC;
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','lsq_blood'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
plot((ODE.x/tRC),gammadot, '.b','LineWidth',2); 



%xticks([0:5:15]);
% yticks([-60]); 
xlim([0 5]) 
xlabel('time', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('gammadot', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{blood}$','$\mu_{blood,0}$', '$\mu_{blood,\infty}$'); 
%set(leg,'Interpreter','latex','Location','southeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.TickLength = [.03 .03];
ax.LineWidth = 1.5;
print('f21New2', '-dpng'); 

%}
