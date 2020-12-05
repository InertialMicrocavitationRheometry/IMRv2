clear; close all; clc;
% reading the data
data = csvread('blood_data.csv',2,0);
skax = (data(data(:,1)~=0,1));
skay = (data(data(:,2)~=0,2));
birox = (data(data(:,3)~=0,3));
biroy = (data(data(:,4)~=0,4));
merix = (data(data(:,5)~=0,5));
meriy = (data(data(:,6)~=0,6));

% Values from Carreau model in Myers et al. PRE
mu0 = 0.056;
mu8 = 0.00345;

% Stacking the three data sets together, can be seperated
gamma = [skax;birox;merix];
mu = [skay;biroy;meriy];
mu = mu(1:end-2);
gamma = gamma(1:end-2);

% defining the carreau equation, x(4) is mu_infinity, x(3) is mu_0, x(2) 
% is n_c, x(1) is \lambda
car = @(x,gamma) x(4) + (x(3)-x(4)).*(1+x(1)^2.*gamma.^2).^((x(2)-1)/2);
% setting initial condition in the guess

% setting lower and upper bounds for the fitting

% running least-squares fit in a nonliner sense, it is still linear
% the arguments in xsol are the values in line 22

%OLD VALUES
%lb = [1,-2,0.1*mu0,0.1*mu8];
%ub = [100,2,2*mu0,2*mu8];
%x0 = [1,0.35,mu0,mu8];
%[xsol,resnorm] = lsqcurvefit(car,x0,gamma,mu,lb,ub);

%NEW VALUES
lb = [1   , 0.2, 0.5*mu0 , 0.5*mu8 ];
ub = [100 , 0.6, 2.0*mu0 , 2*mu8   ];
x0 = (lb+ub)*0.5;
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-8,'FunctionTolerance',1E-8);
[xsol,resnorm] = lsqcurvefit(car,x0,gamma,mu,lb,ub,options);

% plotting the comparisons in a log-log plot to see the comparison more
% fairly, add this figure to the manuscript (after making it prettier)
xl = linspace(0.1,70,100);
figure(1)

loglog(xl,car(xsol,xl), 'k','LineWidth',2);
hold on
loglog(skax,skay,'or','markerfacecolor','r','LineWidth',2);

loglog(birox,biroy,'sb','markerfacecolor','b','LineWidth',2)
loglog(merix,meriy,'g^','markerfacecolor','g','LineWidth',2)

%axis([10E-1 10E1 10E-2 10E-1])
xlim([10E-2 10E1]);
ylim([10E-4 10E-2]);
xticks([10E-2 10E-1 10 10E1]);
yticks([10E-4 10E-3 10E-2]);
xlabel('$\dot{\gamma} = \dot{R}/R$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20); 
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
print('carreauFit1', '-dpng'); 


% plotting the data alongside the fit in a regular figure
figure(2)
hold on;
plot(xl,car(xsol,xl),'-k','LineWidth',2)
hold on
plot(skax,skay,'or','markerfacecolor','r','LineWidth',2);
plot(birox,biroy,'sb','markerfacecolor','b','LineWidth',2)
plot(merix,meriy,'g^','markerfacecolor','g','LineWidth',2)

ylim([0 0.07])
xlabel('$\dot{\gamma} = \dot{R}/R$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20);    
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
print('carreauFit2', '-dpng'); 

%%
%Data sets separated

%SALAK


car_s = @(x,skax) x(4) + (x(3)-x(4)).*(1+x(1)^2.*skax.^2).^((x(2)-1)/2);
lb_s = [1,0.2,0.5*mu0,0.5*mu8];
ub_s = [100,0.6,2*mu0,2*mu8];
x0_s = (lb+ub)*0.5;
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-8,'FunctionTolerance',1E-8);
[xsol_s,resnorm_s] = lsqcurvefit(car_s,x0_s,skax,skay,lb_s,ub_s,options);



xl = linspace(0.1,70,100);
figure(3)
loglog(xl,car_s(xsol_s,xl),skax,skay,'o','LineWidth',2)


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
print('carreauFit1_salak', '-dpng'); 


% plotting the data alongside the fit in a regular figure
figure(4)
hold on;
plot(xl,car_s(xsol_s,xl),'-k','LineWidth',2)
plot(skax,skay,'ob','LineWidth',2)

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
print('carreauFit2_salak', '-dpng'); 



% BIRO

car_b = @(x,birox) x(4) + (x(3)-x(4)).*(1+x(1)^2.*birox.^2).^((x(2)-1)/2);
lb_b = [1,0.2,0.5*mu0,0.5*mu8];
ub_b = [100,0.6,2*mu0,2*mu8];
x0_b = (lb+ub)*0.5;
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-8,'FunctionTolerance',1E-8);
[xsol_b,resnorm_b] = lsqcurvefit(car_b,x0_b,birox,biroy,lb_b,ub_b,options);

xl = linspace(0.1,70,100);
figure(5)
loglog(xl,car_b(xsol_b,xl),birox,biroy,'o','LineWidth',2)


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
print('carreauFit1_biro', '-dpng'); 


% plotting the data alongside the fit in a regular figure
figure(6)
hold on;
plot(xl,car_b(xsol_b,xl),'-k','LineWidth',2)
plot(birox,biroy,'ob','LineWidth',2)

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
print('carreauFit2_biro', '-dpng'); 






% MERRILL

car_m = @(x,merix) x(4) + (x(3)-x(4)).*(1+x(1)^2.*merix.^2).^((x(2)-1)/2);
lb_m = [1,0.2,0.5*mu0,0.5*mu8];
ub_m = [100,0.6,2*mu0,2*mu8];
x0_m = (lb+ub)*0.5;
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-8,'FunctionTolerance',1E-8);
[xsol_m,resnorm_m] = lsqcurvefit(car_m,x0_m,merix,meriy,lb_m,ub_m,options);

xl = linspace(0.1,70,100);
figure(7)
loglog(xl,car_m(xsol_m,xl),merix,meriy,'o','LineWidth',2)


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
print('carreauFit1_meri', '-dpng'); 


% plotting the data alongside the fit in a regular figure
figure(8)
hold on;
plot(xl,car_m(xsol_m,xl),'-k','LineWidth',2)
plot(merix,meriy,'ob','LineWidth',2)

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
print('carreauFit2_meri', '-dpng'); 

%%
figure(9)

loglog(xl,car(xsol,xl), 'k','LineWidth',2);
hold on
loglog(skax,skay,'or','markerfacecolor','r','LineWidth',2);
loglog(xl,car_s(xsol_s,xl),'--r','LineWidth',2)

loglog(birox,biroy,'sb','markerfacecolor','b','LineWidth',2)
loglog(xl,car_b(xsol_b,xl),':b','LineWidth',2)

loglog(merix,meriy,'g^','markerfacecolor','g','LineWidth',2)
loglog(xl,car_m(xsol_m,xl),'-.g','LineWidth',2)

%axis([10E-1 10E1 10E-2 10E-1])
xlim([10E-2 10E1]);
ylim([10E-4 10E-2]);
xticks([10E-2 10E-1 10 10E1]);
yticks([10E-4 10E-3 10E-2]);
xlabel('$\dot{\gamma} = \dot{R}/R$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20); 
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
print('carreauFit3', '-dpng'); 







