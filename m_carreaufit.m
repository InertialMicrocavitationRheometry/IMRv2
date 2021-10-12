clear; close all; clc;
% reading the data
routines = strcat(pwd,'/routines');
addpath(routines);
data = csvread('./blooddata/blood_data.csv',2,0);
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
mu = mu(1:end);
% gamma = gamma(1:end);

% defining the carreau equation, x(4) is mu_infinity, x(3) is mu_0, x(2) 
% is n_c, x(1) is \lambda
% car = @(x,gamma) x(4) + (x(3)-x(4)).*(1+x(1)^2.*gamma.^2).^((x(2)-1)/2);
car = @(x,gamma) x(4) + (x(3)-x(4)).*(1./(1+(x(1).*gamma).^x(2)));
% setting initial condition in the guess
% setting lower and upper bounds for the fitting
lb = [0.2    ,  0.1, mu0 , mu8 ];
ub = [10     ,  2, mu0 , mu8 ];
x0 = (lb+ub)*0.5;
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-8,'FunctionTolerance',1E-8);
% running least-squares fit in a nonliner sense, it is still linear
% the arguments in xsol are the values in line 22
[xsol,resnorm] = lsqcurvefit(car,x0,gamma,mu,lb,ub,options)

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
xlabel('$\dot{\varsigma}$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20); 
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
print('./fitfigs/carreauFit1', '-dpng'); 


% plotting the data alongside the fit in a regular figure
figure(2)
hold on;
plot(xl,car(xsol,xl),'-k','LineWidth',2)
hold on
plot(skax,skay,'or','markerfacecolor','r','LineWidth',2);
plot(birox,biroy,'sb','markerfacecolor','b','LineWidth',2)
plot(merix,meriy,'g^','markerfacecolor','g','LineWidth',2)

ylim([0 0.07])
xlabel('$\dot{\varsigma}$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20);    
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
print('./fitfigs/carreauFit2', '-dpng'); 


% Skalak
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-8);
[xsol_s,resnorm_s] = lsqcurvefit(car,x0,skax,skay,lb,ub,options)

% Biro
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-8);
[xsol_b,resnorm_b] = lsqcurvefit(car,x0,birox,biroy,lb,ub,options)

% Merrill
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-8);
[xsol_m,resnorm_m] = lsqcurvefit(car,x0,merix,meriy,lb,ub,options)

figure(3)
loglog(xl,car(xsol,xl), 'k','LineWidth',2);
hold on
loglog(skax,skay,'or','markerfacecolor','r','LineWidth',2);
loglog(xl,car(xsol_s,xl),'--r','LineWidth',2)

loglog(birox,biroy,'sb','markerfacecolor','b','LineWidth',2)
loglog(xl,car(xsol_b,xl),':b','LineWidth',2)

loglog(merix,meriy,'g^','markerfacecolor','g','LineWidth',2)
loglog(xl,car(xsol_m,xl),'-.g','LineWidth',2)

%axis([10E-1 10E1 10E-2 10E-1])
xlim([10E-2 10E1]);
ylim([10E-4 10E-2]);
xticks([10E-2 10E-1 10 10E1]);
yticks([10E-4 10E-3 10E-2]);
xlabel('$\dot{\varsigma}$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20); 
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
print('./fitfigs/carreauFit3', '-dpng'); 
rmpath(routines);