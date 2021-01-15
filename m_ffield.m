lR = length(Rout(:,1));
lr = rstarres*rstarlim;
R = Rout(:,1);
Rdot = Rout(:,2);
Rdist = Ro_w*ones(lR,lr).*linspace(0,rstarlim,lr);
% calculating the shear
sr = f_sfunc(Rdist,R,Rdot,lR,lr);
mugrid = f_carreau(vmaterial,sr);
munon = (mugrid-mu_inf)/(mu_o-mu_inf);
if mufilter == 1
    [munon,mugrid,sr] = f_mufilter(munon,mugrid,sr,lR,lr);
end
taugrid = -2*mugrid.*sr;
tcon = Tout(1:lR)/tRC;
rcon = Rdist(1,1:lr)'/Ro_w;
[xcon,ycon] = meshgrid(tcon,rcon);
clevels = 50;

% viscosity as a function of r (radial coordinate) plot
figure(6+contourshift)
contourf(xcon',ycon',munon,clevels,'edgecolor','none')
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\mu$*';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
% cbar.Label.Position = [pos(1)*2 pos(2)]; 
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
yticks([0:rstarlimtick:rstarlim])
box on;

% shear rate as a function of r (radial coordinate) plot
figure(7+contourshift)
% plot(tspan,Roft, lm,'LineWidth',2); 
srplot = sr/fnatural;
contourf(xcon',ycon',srplot,clevels,'edgecolor','none')
srmax = max(max(srplot));
srmin = min(min(srplot));
colormap jet;
cbar = colorbar;
set(cbar,'TickLabelInterpreter','latex');
cbar.Label.String = '$\dot\varsigma^*$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
% cbar.Label.Position = [pos(1)*2 pos(2)]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([srmin srmax]);
xlabel('\it{t}$^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{r}$^*$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks([1:1:tmag])
yticks([0:rstarlimtick:rstarlim])
box on;

% viscous shear stress as a function of r (radial coordinate) plot
figure(8+contourshift)
% plot(tspan,Roft, lm,'LineWidth',2); 
tauplot = taugrid/stress;
contourf(xcon',ycon',tauplot,clevels,'edgecolor','none')
taumax = max(max(tauplot));
taumin = min(min(tauplot));
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}^*$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
% cbar.Label.Position = [pos(1)*2 pos(2)]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([taumin taumax]);
xlabel('\it{t}$^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{r}$^*$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks([1:1:tmag])
yticks([0:rstarlimtick:rstarlim])
box on;