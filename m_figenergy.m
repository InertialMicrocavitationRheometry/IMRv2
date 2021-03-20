[E0,LKE,LPE,BIE,VE,TE] = f_energy(Tout,Rout,Pinf,Sd,rho,pGo,kappa,po,vmaterial);
figure()
hold on;
tstar = Tout*fnatural;
plot(tstar,LKE/E0,'--r','LineWidth',2)
plot(tstar,LPE/E0,'-.b','LineWidth',2)
plot(tstar,BIE/E0,':g','LineWidth',2)
plot(tstar,VE/E0,':m','LineWidth',2)
plot(tstar,TE/E0,'k','LineWidth',2)
xlabel('\it{t}$^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{E}^*$', 'Interpreter', 'Latex', 'FontSize', 20); 
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
ylim([0 1.1])
box on;
saveas(gcf,'./analysisfigs/f_energybwall','png')