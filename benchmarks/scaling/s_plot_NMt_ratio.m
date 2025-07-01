% file s_resolution_check.m
% brief contains script to check resolution capability of IMR

% brief This script solves IMR for different spatial resolutions
clc;
clear;
close;

% define parameter vectors and limits
Nt_vec = 2.^(5:8);
Mt_vec = 4*Nt_vec;

% resolved solution
filename = strcat('fd_Nt_',num2str(1024),'_Mt_',num2str(1024));
load(filename);
Rexact = Rvec;

% calculate total number of combinations
total_comb = numel(Nt_vec);
normvec = zeros(total_comb,1);
relnormvec = normvec;

for idx = 1:total_comb
    Nt = Nt_vec(idx);
    Mt = Mt_vec(idx);
    filename = strcat('fd_Nt_',num2str(Nt),'_Mt_',num2str(Mt));
    load(filename);
    normvec(idx) = norm(Rexact-Rvec,2);
    relnormvec(idx) = norm(abs((Rexact-Rvec)./Rexact),2);
end

C = 10 * 120;
D = 2 * (50^2);
y1 = @(x) C ./ x;
y2 = @(x) D ./ (x.^2);
xVals = logspace(0, 4, 200);
yVals1 = y1(xVals);
yVals2 = y2(xVals);

figure(1)
hold on;
plot(Nt_vec,relnormvec,'rs','MarkerFaceColor','r','MarkerSize',10);
loglog(xVals, yVals1,'k-','LineWidth',2);
loglog(xVals, yVals2,'b--','LineWidth',2);

xlim([10^1 10^4])
ylim([10^-1 10^2])
set(gca,'XScale','log')
set(gca,'yScale','log')
box on;
xlabel('$N_t$, $M_t = 4 * N_t$','Interpreter','Latex','FontSize',12);
ylabel('$L_2(|R/R_{\mathrm{exact}} - 1|)$','Interpreter','Latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16);
set(gcf,'color','w');
saveas(gcf,'resolution_error.png')
