fig_tend = 15;
tickrange= [0:3:fig_tend];
lma = 'r';
% SETTING UP THE FIGURES
figure(1)  
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{R} / R_o$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
box on;
plot(t,R, lm,'LineWidth',2); 
saveas(gcf,'./figs/HN/R_T','png')

figure(2)
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\dot{\varsigma}|_R R_o \sqrt{\rho_o/p_{\infty}}$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
% ylim([-1.5 1.5]*1E-5)
box on;
gammadot_R = -2*U./R;
plot(t,gammadot_R,lm,'LineWidth',2); 
saveas(gcf,'./figs/HN/varsigma_T','png')

figure(3)
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$m|_R$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
ylim([0 1.1])
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
ylim([-.1 1.1])
box on;
f = sf_carreau(v_nc,v_lambda,gammadot_R);
% fave = cumsum(f)./[1:length(f)]';
plot(t,f,lm,'LineWidth',2);
saveas(gcf,'./figs/HN/f_T','png')

figure(4)
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}|_R / p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
box on;
tau = 2*(f./DRe+1./Re8).*gammadot_R;
plot(t,tau,lm,'LineWidth',2);
saveas(gcf,'./figs/HN/tau_T','png')