% display result
if plotresult == 1
    subplot(3,1,1);
    plot(t,R,'k','LineWidth',2);
    hold('on');
    ylabel('$R$','Interpreter','Latex','FontSize',12);
    box on;
    if bubtherm == 1
        if medtherm == 0
            b = zeros(size(a));
        end
        subplot(3,1,2);
        hold on;
        box on;
        semilogy(t,abs(a(end,:)),'k-','LineWidth',2);
        semilogy(t,abs(b(end,:)),'b-','LineWidth',2);
        ylabel('$a_N$, $b_M$, $e_N$','Interpreter','Latex','FontSize',12);
        axis([0 t(end) 1e-20 1]);
        set(gca, 'YScale', 'log');
        if masstrans == 0
            leg1 = legend('$a_N$','$b_M$','Location','NorthEast','FontSize',12);
        else
            leg1 = legend('$a_N$','$b_M$','$e_N$','Location','NorthEast','FontSize',12);
        end
        set(leg1,'Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex','FontSize',16);
        set(gcf,'color','w');
    end
    if spectral == 1
        subplot(3,1,3);
        hold on;
        box on;
        semilogy(t,abs(c(end,:)),'k-','LineWidth',2);
        semilogy(t,abs(d(end,:)),'b-','LineWidth',2);
        xlabel('$t$','Interpreter','Latex','FontSize',12);
        ylabel('$c_P$, $d_P$','Interpreter','Latex','FontSize',12);
        set(gca, 'YScale', 'log');
        axis([0 t(end) 1e-20 1]);
        leg1 = legend('$c_P$','$d_P$','Location','NorthEast','FontSize',12);
        set(leg1,'Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex','FontSize',16);
        set(gcf,'color','w');
    end
end
