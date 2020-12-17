lR = length(Rout(:,1));
lr = 500;
R = Rout(:,1);
Rdot = Rout(:,2);
Rdist = Ro_w*ones(lR,lr).*linspace(0,rstarlim,lr);
sr = sfunc(Rdist,R,Rdot,lR,lr);
mugrid = carreau(vmaterial,sr);
taugrid = -2*mugrid.*sr;
tcon = Tout(1:lR)/tRC;
rcon = Rdist(1,1:lr)'/Ro_w;
[xcon,ycon] = meshgrid(tcon,rcon);
clevels = 50;

figure(6)
contourf(xcon',ycon',mugrid/mu_o,clevels,'edgecolor','none')
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\mu$*';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1)*2 -0.066045549576697]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([0 1]);
xlabel('\it{t}$^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{r}$^*$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks([1:1:tmag])
box on;

figure(7)
% plot(tspan,Roft, lm,'LineWidth',2); 
contourf(xcon',ycon',sr,clevels,'edgecolor','none')
colormap jet;
cbar = colorbar;
set(cbar,'TickLabelInterpreter','latex');
cbar.Label.String = '$\dot{\gamma}~[1/s]$';
pos = get(cbar,'Position');
cbar.Label.Position = [1.661309480667114 -10.730386636465376]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([-10 10]);
xlabel('\it{t}$^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{r}$^*$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks([1:1:tmag])
box on;

figure(8)
% plot(tspan,Roft, lm,'LineWidth',2); 
contourf(xcon',ycon',taugrid,clevels,'edgecolor','none')
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}$ [Pa]';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [1.9735118831907 -0.22]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([-0.2 0.2]);
xlabel('\it{t}$^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{r}$^*$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks([1:1:tmag])
box on;