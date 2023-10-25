fig_tend = 30;
tickrange= [0:5:fig_tend];
lR = length(R);
lr_max = 100;
lr_N = 200;
lr_length = 4;
r_coord = ones(lR,lr_N).*logspace(-0.5,lr_length,lr_N);
% calculating the shear
varsigmadot_r = f_gammadot_r(r_coord,R,U,lR,lr_N);
f_r = sf_carreau(v_nc,v_lambda,varsigmadot_r);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
ycon = log10(ycon);
clevels = 50;
feps_r = f_f_filter(f_r,lR,lr_N);
tau_r = 2*(feps_r./DRe+1./Re8).*varsigmadot_r;
min_tau_r = min(min(tau_r))
max_tau_r = max(max(tau_r))
ntau_r = tau_r/max_tau_r;

% f contour figure
figure(5)
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$m(\dot{\varsigma})$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -0.04]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([0 1]);
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
box on;
plot(t, log10(R),lm,'LineWidth',3); 
contourf(xcon',ycon',f_r,clevels,'edgecolor','none')
saveas(gcf,'./figs/baseline/fcon_T','png')

% shear stress contour figure
figure(6)
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/\mathrm{max}(\tau_{rr})$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [1.5*pos(1) -1.2];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([-1 1])
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
box on;
plot(t, log10(R),lm,'LineWidth',3); 
contourf(xcon',ycon',ntau_r,clevels,'edgecolor','none');
saveas(gcf,'./figs/baseline/taucon_T','png')