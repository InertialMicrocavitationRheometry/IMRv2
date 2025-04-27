fig_tend = 30;
tickrange= 0:5:fig_tend;
lR = length(R);
lr_max = 100;
lr_N = 200;
lr_length = 4;
r_coord = ones(lR,lr_N).*logspace(-0.5,lr_length,lr_N);
% calculating the shear
varsigmadot_r = f_gammadot_r(r_coord,R,U,lR,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
ycon = log10(ycon);
clevels = 50;
feps_r = f_f_filter(f_r,lR,lr_N);
tau_r = 2*(feps_r./DRe+1./Re8).*varsigmadot_r;
ntau_r = tau_r/max_tau_r;

% f contour figure
figure(1)
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
clim([0 1]);
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
contourf(xcon',ycon',ntau_r,clevels,'edgecolor','none')
saveas(gcf,'./figs/baseline/fcon_T','png')

% shear as a function of r (radial coordinate) calculation
function sofr = f_gammadot_r(r,R,Rdot,N,M)
    
    sofr = zeros(N,M);
    for i = 1:N
        for j = 1:M
            if r(i,j) >=  R(i,1)
                sofr(i,j) = -2*Rdot(i,1)*(R(i,1)^2)/((r(i,j))^3);
            else
                sofr(i,j) = NaN;
            end
        end
    end
end

function [feps_r] = f_f_filter(f_r,N,M)
    eps = 0.01;
    feps_r = zeros(size(f_r));
    for i = 1:N
        for j = 1:M
            if f_r(i,j) > 1-eps
                feps_r(i,j) = NaN;
            elseif f_r(i,j) < eps
                feps_r(i,j) = NaN;
            else
                feps_r(i,j) = f_r(i,j);
            end
        end
    end
    
end