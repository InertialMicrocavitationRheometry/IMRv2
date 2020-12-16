lR = length(Rout(:,1));
lr = 500;
R = Rout(:,1);
Rdot = Rout(:,2);
Rdist = Ro_w*ones(lR,lr).*linspace(0,8,lr);
sr = sfunc(Rdist,R,Rdot,lR,lr);
mugrid = carreau(vmaterial,sr);
% taugrid = -2*mugrid.*sr;
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
cbar.Label.Position = [pos(1)*2 pos(2)*(-0.3)]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([0 1]);
xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('\it{r*}','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;

% figure(7)
% contourf(xcon',ycon',taugrid,clevels,'edgecolor','none')
% colormap jet;
% cbar = colorbar;
% cbar.Label.String = '$\tau_{rr}$ [Pa]';
% set(cbar,'TickLabelInterpreter','latex');
% pos = get(cbar,'Position');
% cbar.Label.Position = [pos(1)*2 pos(2)*(-3)]; 
% cbar.Label.Rotation = 0;
% cbar.Label.Interpreter = 'latex';
% % caxis([0 1])
% xlabel('\it{t*}', 'Interpreter', 'Latex', 'FontSize', 20); 
% ylabel('\it{r*}','Interpreter','Latex','FontSize',24);
% set(gcf,'color','w');
% set(gca,'FontName','Times','FontSize',20);
% set(gca,'TickLabelInterpreter','latex')
% xa = gca;
% xa.TickLength = [.03 .03];
% xa.LineWidth = 1.5;
% box on;