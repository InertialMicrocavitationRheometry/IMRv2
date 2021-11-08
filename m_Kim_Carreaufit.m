% clearing the workspace, closing figures, and clearing command window
clear; close all; clc;
savedir = strcat(pwd,'/figs/kim_CarreauFit/');
data = csvread('./data/kim_blood_data.csv',2,0);
hemo55x = (data(data(:,1)~=0,1));
hemo55y = (data(data(:,2)~=0,2))/1000;
hemo45x = (data(data(:,3)~=0,3));
hemo45y = (data(data(:,4)~=0,4))/1000;
hemo35x = (data(data(:,5)~=0,5));
hemo35y = (data(data(:,6)~=0,6))/1000;

% Values from Carreau model in Myers et al. PRE
mu0 = 0.056;
mu8 = 0.00345;

% Stacking the three data sets together, can be seperated
gamma = [hemo55x;hemo45x;hemo35x];
mu = [hemo55y;hemo45y;hemo35y];
mu = mu(1:end);
gamma = gamma(1:end);

% defining the carreau equation, x(4) is mu_infinity, x(3) is mu_0, x(2) 
% is n_c, x(1) is \lambda
car = @(x,gamma) mu8 + (mu0-mu8).*(1+x(1)^2.*gamma.^2).^((x(2)-1)/2);
% setting initial condition in the guess
% setting lower and upper bounds for the fitting
lb = [1E-6    ,  0.1 ];
ub = [300   ,   0.8];
x0 = (lb+ub)*0.5;
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-10);
% running least-squares fit in a nonliner sense, it is still linear
% the arguments in xsol are the values in line 22
[xsol,resnorm] = lsqcurvefit(car,x0,gamma,mu,lb,ub,options)

% plotting the comparisons in a log-log plot to see the comparison more
% fairly, add this figure to the manuscript (after making it prettier)
xl = linspace(1,2000,100);

%%% PLOTTING THE X-Y figure
figure(1)
hold on;
plot(xl,car(xsol,xl),'-k','LineWidth',2)
hold on
plot(hemo55x,hemo55y,'or','markerfacecolor','r','LineWidth',2);
plot(hemo45x,hemo45y,'sb','markerfacecolor','b','LineWidth',2)
plot(hemo35x,hemo35y,'g^','markerfacecolor','g','LineWidth',2)
xlabel('$\dot{\varsigma}$ [1/s]','Interpreter','Latex','FontSize',20); 
ylabel('$\mu$ [Pa$\cdot$s]', 'Interpreter', 'Latex', 'FontSize', 20);    
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize', 20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
print(strcat(savedir,'kimFit_xy'), '-dpng'); 

%%% PLOTTING THE LOG-LOG figure
% hemo55
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-10);
[xsol_a,resnorm_s] = lsqcurvefit(car,x0,hemo55x,hemo55y,lb,ub,options)
% hemo45
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-10);
[xsol_b,resnorm_b] = lsqcurvefit(car,x0,hemo45x,hemo45y,lb,ub,options)

% hemo35
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1E-10);
[xsol_c,resnorm_m] = lsqcurvefit(car,x0,hemo35x,hemo35y,lb,ub,options)

figure(2)
loglog(xl,car(xsol,xl), 'k','LineWidth',2);
hold on
loglog(hemo55x,hemo55y,'or','markerfacecolor','r','LineWidth',2);
loglog(xl,car(xsol_a,xl),'--r','LineWidth',2)
loglog(hemo45x,hemo45y,'sb','markerfacecolor','b','LineWidth',2)
loglog(xl,car(xsol_b,xl),':b','LineWidth',2)
loglog(hemo35x,hemo35y,'g^','markerfacecolor','g','LineWidth',2)
loglog(xl,car(xsol_c,xl),'-.g','LineWidth',2)
%axis([10E-1 10E1 10E-2 10E-1])
xlim([5 2E3]);
ylim([1E-3 2E-2]);
xticks([10E-2 10E-1 10 10E1 1E3]);
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
print(strcat(savedir,'kimFit_loglog'), '-dpng'); 