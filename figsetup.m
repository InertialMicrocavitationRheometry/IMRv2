tickrange= [0:tmag];
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
xticks(tickrange)
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
xticks(tickrange)
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
xticks(tickrange)
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;

figure(4)
hold on
xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}$ [Pa]', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
box on;
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xticks(tickrange)
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;

figure(5)
hold on
xlabel('$\dot{\gamma}$ [1/s]', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}$ [Pa]', 'Interpreter', 'Latex', 'FontSize', 20); 
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