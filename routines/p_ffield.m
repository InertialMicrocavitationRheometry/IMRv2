lR = length(Rout(:,1));
lr = rstarres*rstarlim;
R = Rout(:,1);
Rdot = Rout(:,2);
Rdist = Ro_w*ones(lR,lr).*linspace(0,rstarlim,lr);
% calculating the shear
sr = f_sfunc(Rdist,R,Rdot,lR,lr);
vmaterial = 'lsq_blood';
[mugrid,mu_inf,mu_o,~,~] = f_carreau(vmaterial,sr);
munon = (mugrid-mu_inf)/(mu_o-mu_inf);
if mufilter == 1
    [munon,mugrid,sr] = f_mufilter(munon,mugrid,sr,lR,lr,mu_inf);
end
taugrid = -2*mugrid.*sr;
tcon = Tout(1:lR)/tRC;
rcon = Rdist(1,1:lr)'/Ro_w;
[xcon,ycon] = meshgrid(tcon,rcon);
ycon = log10(ycon);
clevels = 50;

% viscosity as a function of r (radial coordinate) plot
figure(6+contourshift)
hold on;
plot(tspan, log10(Roft), lm,'LineWidth',2); 
contourf(xcon',ycon',munon,clevels,'edgecolor','none')
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\mu$*';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -.1]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([0 1]);
xlabel('\it{t}/$t_n$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
xticks([0:5:tmag])
yticks([0:1:log10(rstarlim)])
box on;

% shear rate as a function of r (radial coordinate) plot
figure(7+contourshift)
hold on;
plot(tspan, log10(Roft), lm,'LineWidth',2); 
srplot = sr/fnatural;
contourf(xcon',ycon',srplot,clevels,'edgecolor','none')
colormap jet;
cbar = colorbar;
set(cbar,'TickLabelInterpreter','latex');
cbar.Label.String = '$\dot{\varsigma}/f_n$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
caxis([-2E-5 2E-5]);
cbar.Label.Rotation = 0;
cbar.Label.Position = [pos(1) -2.5E-5]; 
cbar.Label.Interpreter = 'latex';
xlabel('\it{t}/$t_n$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
xticks([0:5:tmag])
% yticks([0:rstarlimtick:rstarlim])
box on;

% viscous shear stress as a function of r (radial coordinate) plot
figure(8+contourshift)
hold on;
tauplot = taugrid/stress;
contourf(xcon',ycon',tauplot,clevels,'edgecolor','none')
plot(tspan,log10(Roft), lm,'LineWidth',2); 
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/B$';
set(cbar,'TickLabelInterpreter','latex');
caxis([-2E-5 2E-5]);
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -2.5E-5]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
xlabel('\it{t}/$t_n$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
xticks([0:5:tmag])
% yticks([0:rstarlimtick:rstarlim])
box on;