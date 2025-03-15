clear;
close all;
clc;
load patterson_waveform_data.mat

data = pattersonwaveformdataall;
clear pattersonwaveformdataall;
data = table2array(data);

%Mostly Linear
t_ML = data(:,1); %time [nondimensional]
p_ML = data(:,2); %input pressure [nondimensional]
rmzeros = (t_ML~=0);
t_ML = t_ML(rmzeros);
p_ML = p_ML(rmzeros);
poly = spline(t_ML,p_ML);
dpoly = fnder(poly,1);
figure(1)
hold on;
% plot(t_ML,p_ML,'.')
trange = linspace(0,15,1E4);
pnew = ppval(poly,trange);
% dpnew = ppval(dp_ML,trange);
plot(trange,pnew,'k','LineWidth',3)
% plot(trange,dpnew)
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$p / p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20);
set(gcf,'color','w'); %Changes background to white
set(gca, 'FontName', 'Times', 'FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
% xticks(tickrange)
% xlim([0 fig_tend])
box on;
% saveas(gcf,'../figs/ML/P_T','png')
save('ml.mat','poly','dpoly');

%Moderately Nonlinear
t_MN = data(:,3); %time [nondimensional]
p_MN = data(:,4); %input pressure [nondimensional]
rmzeros = (t_MN~=0);
t_MN = t_MN(rmzeros);
p_MN = p_MN(rmzeros);
poly = spline(t_MN,p_MN);
dpoly = fnder(poly,1);
figure(2)
hold on;
% plot(t_MN,p_MN,'.');
trange = linspace(0,14,1E5);
pnew = ppval(poly,trange);
% dpnew = ppval(dp_MN,trange);
plot(trange,pnew,'k','LineWidth',3);
% plot(trange,dpnew);
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$p / p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20);
set(gcf,'color','w'); %Changes background to white
set(gca, 'FontName', 'Times', 'FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
% xticks(tickrange)
% xlim([0 fig_tend])
box on;
% saveas(gcf,'../figs/MN/P_T','png')
save('mn.mat','poly','dpoly');

%Highly Nonlinear
t_HN = data(:,5); %time [nondimensional]
p_HN = data(:,6); %input pressure [nondimensional]
rmzeros = (t_HN~=0);
t_HN = t_HN(rmzeros);
p_HN = p_HN(rmzeros);
poly = spline(t_HN,p_HN);
dpoly = fnder(poly,1);
figure(3)
hold on;
% plot(t_HN,p_HN,'.');
trange = linspace(0,14,1E5);
pnew = ppval(poly,trange);
% dpnew = ppval(dp_HN,trange);
plot(trange,pnew,'k','LineWidth',3);
% plot(trange,dpnew);
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$p / p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20);
set(gcf,'color','w'); %Changes background to white
set(gca, 'FontName', 'Times', 'FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
% xticks(tickrange)
% xlim([0 fig_tend])
box on;
% saveas(gcf,'../figs/HN/P_T','png')
save('hn.mat','poly','dpoly');