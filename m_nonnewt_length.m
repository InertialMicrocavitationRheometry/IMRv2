clear all; close all; clc;

dbar = @(muu,mul,nc) (1./(muu.^(2./(nc-1)-1))).^(1/6) - (1./(mul.^(2./(nc-1))-1)).^(1/6);

figure(1)
hold on;

eps = 1E-1;
muu = 1-eps;
mul = eps;
nc = linspace(0,1-1E-6,300);
plot(nc,dbar(muu,mul,nc),'-k','LineWidth',2)

eps = 1E-2;
muu = 1-eps;
mul = eps;
nc = linspace(0,1-1E-6,300);
plot(nc,dbar(muu,mul,nc),'--r','LineWidth',2)

xlabel('\it{$n$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$(r_2-r_1)/\ell$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); % Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.02 .02];
xa.LineWidth = 1.5;
tickrange= [0:0.2:1];
xticks(tickrange)
xlim([0 1.1])
ylim([0 1.1])
box on;
saveas(gcf,'./figs/nonnewt_length/f_nonnewt_length','png')