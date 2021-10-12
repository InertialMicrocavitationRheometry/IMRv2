%%
close;
figure(1)
hold on 
%Carreau model for blood
mu_o = 0.056; 
mu_inf = 0.0345; 
gammadot = linspace(0,5,1E3);
tau_wo = mu_o*gammadot(gammadot<1);
tau_winf = mu_inf*gammadot(gammadot>1.5)+mu_o*0.32;
tau_blood = carreau('blood',gammadot).*gammadot;
plot(gammadot,tau_blood*1e3,'k','LineWidth',3); 
plot(gammadot(gammadot<1),tau_wo*1e3,'b-.','LineWidth',3); 
plot(gammadot(gammadot>1.5),tau_winf*1e3,'r--','LineWidth',3); 
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
%print('f1New2', '-jpg'); 
saveas(gcf,'f1New2.jpg')

%%
close;
figure(2)
hold on
tRC = 0.915*sqrt(rho*Ro_w^2/(A-po-pv));   % Taken from Brennen equation 2.40
%Water
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','water'); 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'b','LineWidth',2); 

%blood inf
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_inf');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'g:','LineWidth',2); 

%Blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); %blood - Carreau 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'r--','LineWidth',2); 
%blood 0
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_knot');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'k-.','LineWidth',2); 

xticks([0:2:20]);
yticks([0:.2:1]); 
axis([0 20 0 1]); 
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
print('f2New2', '-dpng'); 
%saveas(gcf,'f2New2.pdf');

%%
close;
figure(3)
hold on
tRC = 0.915*sqrt(rho*Ro_w^2/(A-po-pv));   % Taken from Brennen equation 2.40
%Water
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','powerlaw','water'); 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'b','LineWidth',2); 

%blood constant value infinity
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','powerlaw','mu_inf');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'g-.','LineWidth',2); 

%Blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','powerlaw','blood'); %blood - Carreau 
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'r--','LineWidth',2);
%blood constant value initial
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','powerlaw','mu_knot');
plot((ODE.x/tRC),(ODE.y(1,:)/Ro_w), 'k-.','LineWidth',2); 

xticks([0:2:20]);
yticks([0:.2:1]);
axis([0 20 0 1]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{R/R$_{0}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$', '$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$'); 
%set(leg,'Interpreter','latex','Location','northeast');
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
print('f3New2', '-dpng'); 

%%
close;
figure(4) 
hold on 
mu_blood = 0.056;
mu_water = 8.9*10E-4; 
mu_o = 0.056; 
%Carreau model for blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); 
plot((ODE.x/tRC), carreau('blood',ODE.y(2,:)./ODE.y(1,:))/mu_blood, 'r','LineWidth',2); 
yline(mu_water/mu_o, 'b','LineWidth',2);
yline(mu_o/mu_o, 'k-.','LineWidth',2)
yline(mu_inf/mu_o, 'g-.','LineWidth',2)
xticks([0:2:20]);
yticks([0:.2:1.2]); 
axis([0 20 0 1.2]); 
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
print('f4new2', '-dpng'); 

%%
close;
figure(5) 
hold on 
%mu_blood = 0.056; 
%Carreau model for blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('blood',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'k','LineWidth',2); 


[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_knot'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('mu_knot',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'--r','LineWidth',2); 

%{
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_inf'); 
gammadot = ODE.y(2,:)./ODE.y(1,:);
tau = carreau('mu_inf',gammadot).*gammadot/1E6;
plot(ODE.x/tRC,tau,'b-.','LineWidth',2); 
%}

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
print('f5New2', '-dpng'); 

%%
close;
figure(6) 
hold on 
mu_blood = 0.056; 
%Carreau model for blood
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','powerlaw','blood');
plot((ODE.x/tRC), powerlaw('blood',ODE.y(2,:)./ODE.y(1,:))/mu_blood, '-r','LineWidth',2);
yline(mu_water/mu_blood, 'b','LineWidth',2);
yline(mu_o/mu_blood, 'k-.','LineWidth',2)
yline(mu_inf/mu_blood, 'g-.','LineWidth',2)
xticks([0:2:20]);
yticks([0:0.01:0.03]); 
axis([0 20 0 .03]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\mu$/$\mu_{0}$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('Power Law'); 
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
print('f6new2', '-dpng'); 

%%
close;
figure(7)
hold on 
%Liver
shearmodulus1 = 1.8e3; %kPa
Ca1 = (rho*c*c)/shearmodulus1;
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca1,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); %blood - Carreau
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'k','LineWidth',2);
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca1,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_knot'); %blood - Carreau
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'k--','LineWidth',2);
%Hepatic Artery
shearmodulus4 = 210e3; %kPa
Ca4 = (rho*c*c)/shearmodulus4;
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca4,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); %
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'b','LineWidth',2);
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca4,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_knot'); %
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'b-.','LineWidth',2);
xticks([0:1:10]);
yticks([0:0.2:1]); 
axis([0 10 0 1]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{R/R$_{0}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('Liver, G = 1.8 kPa, $\mu_{blood}$','Liver, G = 1.8 kPa, $\mu_{blood,0}$',...
%    'Artery, G = 210 kPa, $\mu_{blood}$','Artery, G = 210 kPa, $\mu_{blood,0}$'); 
%set(leg,'Interpreter','latex','Location','northeast');
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
print('f7new2', '-dpng'); 

%%
close; %New for Figure 5
figure(8) 
hold on 
%Gall Bladder
shearmodulus1 = 85e3; %kPa
Ca1 = (rho*c*c)/shearmodulus1;
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca1,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); %blood - Carreau
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'k','LineWidth',2);
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca1,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_knot'); %blood - Carreau
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'k--','LineWidth',2);
%Stomach
shearmodulus4 = 0.637e3; %kPa
Ca4 = (rho*c*c)/shearmodulus4;
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca4,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood'); %
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'b','LineWidth',2);
[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca4,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','mu_knot'); %
plot((ODE.x/tRC),ODE.y(1,:)/Ro_w, 'b-.','LineWidth',2);
xticks([0:1:10]);
yticks([0:0.2:1]); 
axis([0 10 0 1]); 
xlabel('\it{t/t$_{RC}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{R/R$_{0}$}', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('Gall Bladder, G = 85 kPa, $\mu_{blood}$','Gall Bladder, G = 85 kPa, $\mu_{blood,0}$',...
%   'Stomach, G = 0.64 kPa, $\mu_{blood}$','Stomach, G = 0.64 kPa, $\mu_{blood,0}$'); 
%set(leg,'Interpreter','latex','Location','northeast');
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
%imwrite(figure(8),'f8New.png')
print('f8New2', '-dpng'); 
%saveas(gcf,'f8new2','epsc');
%%
close; 
figure(9) 
hold on 

data = readmatrix('blood_data.csv', 'HeaderLines', 2);

salak_X = data(1:8,1);
salak_Y = data(1:8,2);
biro_X = data(:,3);
biro_Y = data(:,4);
merrill_X = data(1:9,5);
merrill_Y = data(1:9,6);


plot(salak_X,salak_Y,'ob','LineWidth',2)
hold on
plot(biro_X,biro_Y,'ok','LineWidth',2)
plot(merrill_X,merrill_Y,'oc','LineWidth',2)


gamma_dot = [salak_X; biro_X; merrill_X];
mu_data = [salak_Y; biro_Y; merrill_Y];
mu_inf = 0.0045; 
muo = 0.056; 
lambda = 3.313; 
%nc = 0.3568; 


% AVG ALL DATA SETS
mu_c = @(x,gamma_dot) mu_inf + (muo - mu_inf).*(1+(lambda).^2*(gamma_dot).^2).^((x(1)-1)./2); 
%epsilion_0 = sqrt(sum(abs(mu_carreau-mu_data).^2)) ./ sqrt(sum(mu_data.^2));

fun = @(x) sqrt(sum(abs(mu_c(x,gamma_dot)-mu_data).^2)) ./ sqrt(sum(mu_data.^2));
x0 = [2];
nc = fminsearch(fun,x0)

mu_carreau = mu_inf + (muo - mu_inf).*(1+(lambda).^2*(gamma_dot).^2).^((nc-1)./2); 
plot(gamma_dot, mu_carreau, 'or','LineWidth',2); 





%axis([0 60 0 0.1])
leg = legend('Salak', 'Biro', 'Merrill', 'Carreau1'); 
set(leg,'Interpreter','latex','Location','northeast');
set(leg,'FontSize',18);
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
print('f9new2', '-dpng'); 

%%
close; 
figure(10) 
hold on 

data = readmatrix('blood_data.csv', 'HeaderLines', 2);

salak_X = data(1:8,1);
salak_Y = data(1:8,2);
biro_X = data(:,3);
biro_Y = data(:,4);
merrill_X = data(1:9,5);
merrill_Y = data(1:9,6);


plot(salak_X,salak_Y,'ob','LineWidth',2)
hold on
plot(biro_X,biro_Y,'ok','LineWidth',2)
plot(merrill_X,merrill_Y,'oc','LineWidth',2)


gamma_dot = [salak_X; biro_X; merrill_X];
mu_data = [salak_Y; biro_Y; merrill_Y];
mu_inf = 0.0045; 
muo = 0.056; 
lambda = 3.313; 
%nc = 0.3568; 


% AVG INDIVIDUALLY
mu_c_S = @(x,salak_X) mu_inf + (muo - mu_inf).*(1+(lambda).^2*(salak_X).^2).^((x(1)-1)./2); 
fun = @(x) sqrt(sum(abs(mu_c(x,salak_X)-salak_Y).^2)) ./ sqrt(sum(salak_Y.^2));
x0 = [2];
nc_S = fminsearch(fun,x0)

mu_c_B = @(x,biro_X) mu_inf + (muo - mu_inf).*(1+(lambda).^2*(biro_X).^2).^((x(1)-1)./2); 
fun = @(x) sqrt(sum(abs(mu_c(x,biro_X)-biro_Y).^2)) ./ sqrt(sum(biro_Y.^2));
x0 = [2];
nc_B = fminsearch(fun,x0)

mu_c_M = @(x,merrill_X) mu_inf + (muo - mu_inf).*(1+(lambda).^2*(merrill_X).^2).^((x(1)-1)./2); 
fun = @(x) sqrt(sum(abs(mu_c(x,merrill_X)-merrill_Y).^2)) ./ sqrt(sum(merrill_Y.^2));
x0 = [2];
nc_M = fminsearch(fun,x0)

nc_avg = mean([nc_S;nc_B;nc_M])
mu_carreau_avg = mu_inf + (muo - mu_inf).*(1+(lambda).^2*(gamma_dot).^2).^((nc_avg-1)./2); 
plot(gamma_dot, mu_carreau_avg, 'og','LineWidth',2); 

%axis([0 60 0 0.1])
leg = legend('Salak', 'Biro', 'Merrill', 'Carreau2'); 
set(leg,'Interpreter','latex','Location','northeast');
set(leg,'FontSize',18);
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
print('f10new2', '-dpng'); 
%%
close; 
figure(11) 
hold on 

data = readmatrix('blood_data.csv', 'HeaderLines', 2);

salak_X = data(1:8,1);
salak_Y = data(1:8,2);
biro_X = data(:,3);
biro_Y = data(:,4);
merrill_X = data(1:9,5);
merrill_Y = data(1:9,6);


plot(salak_X,salak_Y,'ob','LineWidth',2)
hold on
plot(biro_X,biro_Y,'ok','LineWidth',2)
plot(merrill_X,merrill_Y,'oc','LineWidth',2)

[ODE]=bubblewall_solver(Ro_w,A,w,We_w,Ca,rho,po,pv,c,tfinal,'KM','NeoH','Carreau','blood');

gammadot = linspace(0,70,length(ODE.x));

%np = 0.6;
%K = 0.035;
K = 0.028;
np = 0.28;

mu_powerlaw = K * gammadot.^(np-1)+1E-6;
plot(gammadot,mu_powerlaw, '-k','LineWidth',2);
%axis([0 60 0 0.1])
leg = legend('Salak', 'Biro', 'Merrill', 'Power Law'); 
set(leg,'Interpreter','latex','Location','northeast');
set(leg,'FontSize',18);
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
print('f11new2', '-dpng'); 