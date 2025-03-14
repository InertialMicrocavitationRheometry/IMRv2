% file s_resolution_check.m
% brief contains script to check resolution capability of IMR

% brief This script solves IMR for different spatial resolutions
clc;
clear;
close;

load('fd2048.mat');
load('fd_ratio_scaling.mat');

% define parameter vectors and limits
Mt_vec = 2.^(5:10);
Nt_vec = 2.^(4:9)+1;

% calculate total number of combinations
total_comb = numel(Nt_vec);
normvec = zeros(total_comb,1);
relnormvec = normvec;
dtvec = [0;diff(tfd2048)];

for idx = 1:total_comb
    normvec(idx) = norm((Rfd2048-Rvec{idx}(:)),2);
    relnormvec(idx) = norm(abs((Rfd2048-Rvec{idx}(:))./Rfd2048),2);
end
Nvec = Nt_vec+0*Mt_vec;


C = 30 * 120;  
D = 1 * (50^2);  
y1 = @(x) C ./ x;
y2 = @(x) D ./ (x.^2);
xVals = logspace(0, 4, 200);  
yVals1 = y1(xVals);
yVals2 = y2(xVals);

figure(1)
hold on;
plot(Nvec,relnormvec,'rs','MarkerFaceColor','r','MarkerSize',7);
% plot(Nvec,normvec,'^k','MarkerFaceColor','k');
loglog(xVals, yVals1,'k-','LineWidth',2); 
loglog(xVals, yVals2,'b--','LineWidth',2);

xlim([10^0 10^4])
ylim([10^-1 10^3])
set(gca,'XScale','log')
set(gca,'yScale','log')
box on;
xlabel('$N_t$, $M_t = 2 * N_t$','Interpreter','Latex','FontSize',12);
ylabel('$L_2(|R/R_e - 1|)$','Interpreter','Latex','FontSize',12);
% leg1 = legend('$\epsilon = |R-R_e|/|R_e|$','$\epsilon = |R-R_e|$','Location','NorthEast','FontSize',12);
% leg1 = legend('$\epsilon = |R/R_e - 1|$','Location','NorthEast','FontSize',12);
% set(leg1,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
set(gcf,'color','w');
saveas(gcf,'resolution_error.png')


% figure(2)
% hold on;
% plot(tfd2048,Rfd2048,'sk');
% legvec = cell(total_comb/2,1);
% legvec{1} = '$N_t = 500$';
% counter = 1;
% for idx = 1:total_comb
%     plot(tvec{idx}(:),Rvec{idx}(:),'--')
%     legvec{counter+1} = strcat('$N_t = ',num2str(Nt_vec(idx)),'$');
%     counter = counter + 1;
% end
% box on;
% xlabel('$t/t_c$','Interpreter','Latex','FontSize',12);
% ylabel('$R(t)/R_o$','Interpreter','Latex','FontSize',12);
% 
% leg1 = legend(legvec,'FontSize',12);
% set(leg1,'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex','FontSize',16);
% set(gcf,'color','w');
% saveas(gcf,'resolution_Roft.png')
