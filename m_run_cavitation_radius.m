% This script can be used as a template to run RP_Cav
clear; close all; clc;
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
R_range = R0*[0.25, 0.5, 1, 2, 4];
p_l_mu_max = zeros(length(R_range));
p_l_mu_avg = zeros(length(R_range));
R_max = zeros(length(R_range));
R_ave = zeros(length(R_range));
tfactor = R_range./R0;

for i=1:length(R_range)
    [Pmt] = f_call_parameters(R_range(i),vmaterial);
    Pinf = Pmt(19); v_lambda = Pmt(27);
    Pext_Amp_Freq = [A, f]; 
    tend = 8E-6/tfactor(i);
    [ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
    (tend,R_range(i),NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
    Dim,comp,vmaterial,vmodel);
    vec = f_char_length(R,U,v_lambda)*R_range(i);
    p_l_mu_max(i) = max(vec);
    p_l_mu_avg(i) = mean(vec); 
    R_max(i) = max(R)*R_range(i);
    R_ave(i) = mean(R)*R_range(i);
end

figure(1)
hold on;
plot(R_range/R0,p_l_mu_max,'or','MarkerFaceColor','r','MarkerSize',8);
plot(R_range/R0,R_max,'ob','MarkerFaceColor','b','MarkerSize',8);
plot(R_range/R0,p_l_mu_avg,'^r','MarkerFaceColor','r','MarkerSize',8);
plot(R_range/R0,R_ave,'^b','MarkerFaceColor','b','MarkerSize',8);
xlabel('$R /R_{o}$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$ \ell $', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gca, 'YScale', 'log','XScale', 'log')
set(gcf,'color','w'); 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex');
xlim([1E-1 1E1])
ylim([1E-8 1E-2])
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
saveas(gcf,'./figs/param/ell_radius','png')

rmpath(routines);
rmpath(post);