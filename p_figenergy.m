[E0,EPD,LKE,LPE,BIE,VE,TE] = f_energy(Tout,Rout,Pinf,DPinf,Sd,rho,pGo,kappa,po,...
    vmaterial,force);
dem = E0;
figure(6)
hold on;
tstar = Tout*fnatural;
plot(tstar,LKE./dem,':r','LineWidth',2)
plot(tstar,LPE./dem,'-.b','LineWidth',2)
plot(tstar,BIE./dem,'--g','LineWidth',2)
plot(tstar,(VE-EPD)./dem,'-om','LineWidth',2,'MarkerIndices',1:ceil(length(tstar)/20):length(tstar));
% plot(tstar,EPD./dem,'-oc','LineWidth',2,'MarkerIndices',1:ceil(length(tstar)/20):length(tstar));
plot(tstar,TE./dem,'k','LineWidth',2)
xlabel('$\it{t}/t_{n}$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{E}/E_o$', 'Interpreter', 'Latex', 'FontSize', 20); 
%leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 tmag])
box on;
fn = strcat('./analysisfigs/f_energybwall',filesuffix);
saveas(gcf,fn,'png')