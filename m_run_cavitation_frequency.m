% This script can be used as a template to run RP_Cav
% clear; close all; clc;
% Adding the routines path
routines = strcat(pwd,'/src');
post = strcat(pwd,'/post_processing');
addpath(genpath(routines));
addpath(post);

% All quantities below are in SI units (s,m)
% time to run simulation % tend = 2E-6;
% Ammount of nodes inside the bubble (~100 is a good to start)
NT = 5; % NT = 50; 
% Ammount of nodes outside the bubble (~100 is a good to start)
NTM = 5; % NTM = 60; 
% Type of external forcing. Refer to RP_Cav
Pext_type = 'sn'; 
% 1 : display simulation time, 0 : do not display
disptime = 1; 
% Thermal effects inside bubble, 1: yes, 0: no
Tgrad = 0; 
% Thermal effects outside bubble, 1: yes, 0: no
Tmgrad = 0;
% Vapor diffusion effects, 1: yes, 0: no
Cgrad = 0;  
% Output variables in dimensional form, 1: yes, 0: no
Dim = 0; 
% Activates the effect of compressibility 
% 1: Keller-Miksis w/ pressure, 0: Rayleigh-Plesset
comp = 1; 
% material to calculate viscosity
vmaterial = 'blood_combined';
% non-Newtonian model for viscosity
vmodel = 'carreau';
% Baseline values
A = 1E6; 
f = 2E6;
R0= 0.5E-6;

% RUNNING FREQUENCY VARIATION
f_range = f*2.^(-2.5:0.5:2.5);
f_l_mu_max = zeros(length(f_range),1);
f_l_mu_avg = zeros(length(f_range),1);
R_max = zeros(length(f_range),1);
R_ave = zeros(length(f_range),1);
tfactor = f_range./f_range(1);

[Pmt] = f_call_parameters(R0,vmaterial);
t0 = Pmt(14); v_lambda = Pmt(27);
Pinf = Pmt(19);
xrange = f_range*t0*(Pinf/A);
% figure(2)
% hold on;
for i=1:length(f_range)
    Pext_Amp_Freq = [A, f_range(i)]; 
    tend = 2E-6/tfactor(i);
    [ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
    (tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
    Dim,comp,vmaterial,vmodel);
    vec = f_char_length(R,U,v_lambda);
    f_l_mu_max(i) = max(vec);
    f_l_mu_avg(i) = mean(vec); 
    R_max(i) = max(R);
    R_ave(i) = mean(R);
%     figure(2);
%     plot(t,R/R_max(i));
end

yy = [10^4 10^2];
xx = [2*10^-3 2*10^-1];

addpath(genpath(routines));
addpath(post);
figure(1)
hold on;
plot(xx,yy,'k','LineWidth',2);
plot(xrange,f_l_mu_max,'sr','MarkerFaceColor','r','MarkerSize',8);
plot(xrange,R_max,'sb','MarkerFaceColor','b','MarkerSize',8);
% plot(f_range*t0,f_l_mu_avg,'^r','MarkerFaceColor','r','MarkerSize',8);
% plot(f_range*t0,R_ave,'^b','MarkerFaceColor','b','MarkerSize',8);
xlabel('$f R_o\sqrt{\rho_{\ell}/p_{\infty}}$', 'Interpreter', 'Latex', 'FontSize', 20); 
% ylabel('$ \ell $', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gca, 'YScale', 'log','XScale', 'log')
set(gcf,'color','w'); 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex');
% ylim([1E-8 1E-2])
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
% saveas(gcf,"./figs/param/ell_frequency",'png')

% rmpath(routines);
% rmpath(post);