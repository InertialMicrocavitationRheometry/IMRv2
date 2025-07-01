% file s_resolution_check.m
% brief contains script to check resolution capability of IMR

% brief This script solves IMR for different spatial resolutions
clc;
clear;
close;

% define parameter vectors and limits
Mt = 1024;
Nt_vec = 2.^(5:8);

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
    filename = strcat('fd_Nt_',num2str(Nt),'_Mt_',num2str(Mt));
    load(filename);
    normvec(idx) = norm(Rexact-Rvec,2);
    relnormvec(idx) = norm((Rexact-Rvec)./Rexact,2);
end

C = 5 * 120;
D = 1 * (4^2);
y1 = @(x) C ./ x;
y2 = @(x) D ./ (x.^2);
xVals = logspace(-3, 6, 200);
yVals1 = y1(xVals);
yVals2 = y2(xVals);

figure(1)
hold on;
plot(Nt_vec,relnormvec,'rs','MarkerFaceColor','r','MarkerSize',10);
plot(Nt_vec,normvec,'^k','MarkerFaceColor','k','MarkerSize',10);
loglog(xVals, yVals1,'k-','LineWidth',2);
loglog(xVals, yVals2,'b--','LineWidth',2);
xlim([10^0 10^4])
ylim([10^-2 10^2])
set(gca,'XScale','log')
set(gca,'yScale','log')
box on;
xlabel('$N_t$','Interpreter','Latex','FontSize',12);
ylabel('$L_2(|R/R_{\mathrm{exact}} - 1|)$','Interpreter','Latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16);
set(gcf,'color','w');
saveas(gcf,'Nt_error.png')

figure(2)
hold on;
plot(Rexact,'sk');
legvec = cell(total_comb,1);
legvec{1} = '$N_t = 1000$';
counter = 1;
for idx = 1:total_comb
    Nt = Nt_vec(idx);
    filename = strcat('fd_Nt_',num2str(Nt),'_Mt_',num2str(Mt));
    load(filename);
    plot(Rvec,'--')
    legvec{counter+1} = strcat('$N_t = ',num2str(Nt_vec(idx)),'$');
    counter = counter + 1;
end
box on;
xlabel('$t/t_c$','Interpreter','Latex','FontSize',12);
ylabel('$R(t)/R_o$','Interpreter','Latex','FontSize',12);

leg1 = legend(legvec,'FontSize',12);
set(leg1,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
set(gcf,'color','w');
