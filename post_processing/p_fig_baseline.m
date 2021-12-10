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
ylabel('$\dot{\varsigma}^|_R$', 'Interpreter', 'Latex', 'FontSize', 20); 
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
box on;
[Pmt] = f_call_parameters(R0,vmaterial);
v_nc = Pmt(26); v_lambda = Pmt(27);
v_mu8 = Pmt(28); v_Dmu = Pmt(29);
f = sf_carreau(v_nc,v_lambda,gammadot_R);
fave = cumsum(f)./[1:length(f)]';
plot(t,f,lm,'LineWidth',2);
plot(t,fave,lma,'LineWidth',2);

figure(4)
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}|_R$', 'Interpreter', 'Latex', 'FontSize', 20); 
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
tau = 2*(f*v_Dmu+v_mu8).*gammadot_R;
plot(t,tau,lm,'LineWidth',2);

figure(5)
hold on
xlabel('$\dot{\varsigma}^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\tau_{rr}|_R$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xlim([0 1])
% ylim([0 1])
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
[sort_gammadot_R,sortI] = sort(gammadot_R);
sort_tau = tau(sortI);
plot(sort_gammadot_R,sort_tau,'LineWidth',2);