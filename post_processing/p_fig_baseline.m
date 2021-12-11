tend = t(end);
tickrange= [0:10:tend];
lm = 'k';
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
xlim([0 tend])
box on;
plot(t,R, lm,'LineWidth',2); 

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
xlim([0 tend])
% ylim([-1.5 1.5]*1E-5)
box on;
gammadot_R = -2*U./R;
plot(t,gammadot_R,lm,'LineWidth',2); 

figure(3)
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{f}|_R$', 'Interpreter', 'Latex', 'FontSize', 20); 
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
xlim([0 tend])
ylim([-.1 1.1])
box on;
f = sf_carreau(v_nc,v_lambda,gammadot_R);
% fave = cumsum(f)./[1:length(f)]';
plot(t,f,lm,'LineWidth',2);
% plot(t,fave,lma,'LineWidth',2);

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
xlim([0 tend])
box on;
tau = 2*(f./DRe+1./Re8).*gammadot_R;
plot(t,tau,lm,'LineWidth',2);

figure(5)
hold on
xlabel('$\dot{\varsigma}|_R R_o \sqrt{\rho_o/p_{\infty}}$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}|_R / p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
axis tight;
xa = gca;
xa.TickLength = [.025 .025];
xa.LineWidth = 1.5;
set(gca,'XLim',[1e-8 1e2],'XTick',10.^(-8:2:2), ...
        'YLim',[1e-10 1e0],'YTick',10.^(-10:2:0))
box on;
sgammadot_R = linspace(1E-8,1E2,100);
sf = sf_carreau(v_nc,v_lambda,sgammadot_R);
tau8 = 2*(0./DRe + 1./Re8)*sgammadot_R;
tau0 = 2*(1./DRe+1./Re8)*sgammadot_R;
tau  = 2*(sf./DRe+1./Re8).*sgammadot_R;
set(gca,'XScale','log','YScale','log')
plot(sgammadot_R,tau,lm,'LineWidth',2);
plot(sgammadot_R,tau8,'b--','LineWidth',2);
plot(sgammadot_R,tau0,'r-.','LineWidth',2);