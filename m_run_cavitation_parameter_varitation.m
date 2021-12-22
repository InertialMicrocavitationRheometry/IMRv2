% This script can be used as a template to run RP_Cav
clear; close all; clc;
% Adding the routines path
routines = strcat(pwd,'\src');
post = strcat(pwd,'\post_processing');
addpath(genpath(routines));
addpath(post);

% All quantities below are in SI units (s,m)
% time to run simulation 
tend = 1.5E-6;
% Initital Radii 
R0 = 0.5E-6;
% Ammount of nodes inside the bubble (~100 is a good to start)
NT = 50; 
% Ammount of nodes outside the bubble (~100 is a good to start)
NTM = 60; 
% Type of external forcing. Refer to RP_Cav
Pext_type = 'sn'; 
% [ Pressure ; Freq ] 
Pext_Amp_Freq =[1E6 2E6]; 
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
%%
% RUNNING BASELINE CODE
[ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27);
lm = 'k';
r_fig_baseline;

% RUNNING BASELINE CODE
vmaterial = 'blood';
[ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
(tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
Dim,comp,vmaterial,vmodel);
[Pmt] = f_call_parameters(R0,vmaterial);
Re8 = Pmt(6); DRe = Pmt(24);
v_nc = Pmt(26); v_lambda = Pmt(27);
lm = 'r--';
r_fig_baseline;
%%
% Baseline values
A = 1E6; 
f = 2E6;
Ro_w= 0.5E-6;

Pext_type = 'sn'; 

% RUNNING FREQUENCY VARIATION
f_range = [.1E6, .5E6, 1E6, 2E6, 5E6];

%f_p_range = ones(1,length(f_range))*A;
%f_Ro_range = ones(1,length(f_range))*Ro_w;
f_l_mu_max = zeros(1,length(f_range));
f_l_mu_avg = zeros(1,length(f_range));



f_names = ['./figs/param/freq/1E5_R_T'; './figs/param/freq/1E5_varsigma_T'; ...
    './figs/param//freq/1E5_f_T'; './figs/param//freq/1E5_tau_T'; ...
    './figs/param/freq/2E5_R_T'; './figs/param/freq/2E5_varsigma_T'; ...
    './figs/param//freq/2E5_f_T'; './figs/param//freq/2E5_tau_T'; ...
    './figs/param/freq/1E6_R_T'; './figs/param/freq/1E6_varsigma_T'; ...
    './figs/param//freq/1E6_f_T'; './figs/param//freq/1E6_tau_T'; ...
    './figs/param/freq/2E6_R_T'; './figs/param/freq/2E6_varsigma_T'; ...
    './figs/param//freq/2E6_f_T'; './figs/param//freq/2E6_tau_T'; ...
    './figs/param/freq/5E6_R_T'; './figs/param/freq/5E6_varsigma_T'; ...
    './figs/param//freq/5E6_f_T'; './figs/param//freq/5E6_tau_T';];

lm_vec = ["r","--b","-.k",":m","--g"]
    
f_name_count = 0;  

vmaterial = 'blood_combined';
for i=1:length(f_range)
    Pext_Amp_Freq = [A, f_range(i)]; 
    [ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
    (tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
    Dim,comp,vmaterial,vmodel);
    [Pmt] = f_call_parameters(R0,vmaterial);
    Re8 = Pmt(6); DRe = Pmt(24);
    v_nc = Pmt(26); v_lambda = Pmt(27);
    lm = lm_vec(i); 
    
    r_fig_variation;
    
    vec = f_char_length(R,U,v_lambda);
    f_l_mu_max(i) = max(vec.');
    f_l_mu_avg(i) = mean(vec.'); 
end

figure()
%semilogx(f_range(3:end),f_l_mu_max(3:end),'or');
semilogx(f_range,f_l_mu_max,'or');
xlabel('$f$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$ \it{l}_{\mu,max} $', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
%xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
saveas(gcf,'./freq/l_mu_max','png')

figure()
%semilogx(f_range(3:end),f_l_mu_max(3:end),'or');
semilogx(f_range,f_l_mu_avg,'or');
xlabel('$f$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$ \it{l}_{\mu,avg} $', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
%xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
saveas(gcf,'./freq/l_mu_avg','png')




