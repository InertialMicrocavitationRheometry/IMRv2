close all; clear; clc;
routines = strcat(pwd,'/src');
addpath(genpath(routines));
R0 = 500E-6; vmaterial = 'blood_combined';
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27);
lm = 'k';

figure(1)
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
set(gca,'XLim',[1e-7 1e0],'XTick',10.^(-6:2:0), ...
        'YLim',[1e-10 1e-2],'YTick',10.^(-10:2:-2))
box on;
sgammadot_R = logspace(-8,0,100);
sf = sf_carreau(v_nc,v_lambda,sgammadot_R);
tau8 = 2*(0./DRe + 1./Re8)*sgammadot_R;
tau0 = 2*(1./DRe+1./Re8)*sgammadot_R;
tau  = 2*(sf./DRe+1./Re8).*sgammadot_R;
set(gca,'XScale','log','YScale','log')
plot(sgammadot_R,tau,lm,'LineWidth',2);
plot(sgammadot_R,tau8,'b--','LineWidth',2);
plot(sgammadot_R,tau0,'r-.','LineWidth',2);

saveas(gcf,'./figs/shear_analysis/tau_shear','png')

figure(2)
hold on;
xlabel('$\dot{\varsigma} R_o \sqrt{\rho_o/p_{\infty}}$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{f}(\dot{\varsigma})$', 'Interpreter', 'Latex', 'FontSize', 20); 
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
set(gca,'XLim',[1e-8 1e2],'XTick',10.^(-8:2:2));
set(gca,'XScale','log');
ylim([-.1 1.1]);
box on;
sgammadot_R = logspace(-8,2,100);
sf = sf_carreau(v_nc,v_lambda,sgammadot_R);
s8 = zeros(size(sf));
s0 = ones(size(sf));
plot(sgammadot_R,sf,lm,'LineWidth',2);
plot(sgammadot_R,s8,'b--','LineWidth',2);
plot(sgammadot_R,s0,'r-.','LineWidth',2);

saveas(gcf,'./figs/shear_analysis/f_shear','png')
rmpath(routines);