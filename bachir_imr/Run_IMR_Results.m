% close all;clc;
warning('off','all')

% Load Data
% load('C:\Users\bachi\Dropbox (University of Michigan)\Bachir\ODV\Data\PA_Calibrated.mat')
% load('/Users/bachirabeid/Dropbox (University of Michigan)/Bachir/ODV/Data/All_Combined.mat')
% load('matlab2.mat')
load('matlab2.mat')

% saveFigs = 'C:\Users\bachi\Dropbox (University of Michigan)\Bachir\ODV\IMR_Modeling\Figures\figs\';
% savePng = 'C:\Users\bachi\Downloads\Temp';%E:\LDV\Temp_figures\';

tau = 0.1E-7;
alp = 0.9;

relax_type = 'none'; % Options: 'none', 'maxwell', 'fd'

% Define RO parameters:
% NOTE: set dimensionless moduli, so that effect is similar at all size
Ca_try = 1E2; 
Re_try = 1E2; % Actually Ze for FD, but size similar

nruns = 3; % 2x3

% Generate save name:

% Define pre-load & spatial mesh (for Maxwell/FD)

npre = 40;
%IX = zeros(npre,4); % MATLAB recommended against pre-allocation ... 
IY = zeros(npre,4); % Only pass time ratio, sin, cos, do rest case-by-case
for kk = 1:(npre)
    % kt = -TX*(1+npre-kk)/npre;
    % xscale = 1.2;
    % kt = (xscale^(1-kk)); 
    kt = (npre+1-kk)/npre;

    % April 2022 version: (L-1) = (Lmax-1)*(t/TX)^3.665
    IY(kk,:) = [kt,kt^3.665,3.665*(kt^2.665),(3.665*2.665)*(kt^1.665)];
end

% Note: Removed t = 0 from IX, since that step is evaluated during ODE solution

% Load meshing
load mesh.mat % Get 'RMesh'

% Other basic parameters
% (We can directly treat them here ...)
P_inf = 101325; % (Pa) Atmospheric Pressure
rho = 1064; % (Kg/m^3) Material Density: PA - 998.2
Uc = sqrt(P_inf/rho); % Characteristic velocity

% Additional parameters needed to find initial gas pressure - Nov.'21
% (!!!--Make sure this is consistent with IMRcall_parameters.m--!!!)
gam_st = 0.072; % (N/m) Water Surface Tension (Called "S" in other function)
T_inf = 298.15; % (K) Far field temp. 

Tgrad = 1; %1 Temperature transfer bt bubble and material
Cgrad = 1;%1;  %1 Vapor-non-condensible gas diffusion
Tmgrad = 1; %0 Off means cold liquid assumption

% Note: P_guess depends on bubble size, so we will treat it inside solver
% loop ...

% Run IMR
FigName ='Collapse_time_PA_F_6um_pfp';
range = [181:191,170:180];%[78:98,99:105];%;
% 
n = 219;
FigName ='Collapse_time_PAA_B_6um';
range = [160:169,148:159];%firbin 0.2% 
n = 169;
% Gelatin Native:                   5% --> 1:6 ; 10% --> 7:15 ; 15% --> 16:22
% Gelatin w/ 12um Droplets:          5% 55:63 ; 10% --> 64:69 ; 15% --> 70:77
% Gelatin 10% w/ 3um & 6um Droplets: 3um --> 40:49 ; 6um --> 50:54
% Gelatin 10% w/ pfh/pfo:            pfh: 23:31 ; pfo: 32:39

% Fribrin Native:           0.2% --> 78:98 ; 1% --> 106:115 ; 4% --> 125:135
% Fribrin w/ 6um Droplets:  0.2% --> 99:105 ; 1% --> 116:124 ; 4% --> 136:147

% PAA Native:           B --> 160:169 ; D --> 181:191 ; E --> 204:213 ; F --> 223:234
% PAA w/ 6um Droplets:  B --> 148:159 ; D --> 170:180 ; E --> 192:203 ; F --> 214:222

% Storage for post-processing
T1X = zeros(nruns,1);

for i =range(1)% 1:length(expts)
    
    G = expts(i).G_best*0.25;
    mu = expts(i).mu_best*1.29 ; %0.2--> 0.85_1.15 ; 1% __> 0.75_1.25
    % expts(i).G_best = 21625;
    % expts(i).mu_best = 0.011;
    eqR = expts(i).Req*expts(i).scale;
    t =expts(i).t;  
    R =expts(i).Roft;
    R0 =expts(i).R0;
    % eqR = R0;
    fileName = expts(i).FileName;

    P_guess = (P_inf+2*gam_st/eqR-Pvsat(T_inf))*(eqR/R0)^(3);
    Lmax = R0/eqR; % Max. stretch

    % Set up simulation duration:
    tspan = R0/3; %/3; % Make a bit longer, for good comparison
    TX = R0/12.6; % Empirically selected duration for initial expansion, matching experiment ...

    % Dimensionalize preload:
    IX = [-IY(:,1)*TX,R0 - IY(:,2)*(R0-eqR),IY(:,3)*(R0-eqR)/TX,-IY(:,4)*(R0-eqR)/(TX^2)];

    
    switch relax_type
        case 'maxwell'
            model = 'zzzen';
            G1 = mu/tau;
        case 'fd'
            model = 'fdkv';
            G1 = alp;
        case 'none'
            model = 'neoHook';
            G1 = inf;
    end


    % Treat bubble content: (Left original code in place)
    NT = 50; % Mesh points in bubble, resolution should be >=500
    NTM = 50; % Mesh points in the material, should be 10
    Pext_type = 'IC'; %'IC' for Flynn, 'ga' for gaussian bubble growth

    if strcmp(Pext_type,'IC')
        Pext_Amp_Freq = [P_guess; 0];%[226; 0]; % Tune first number to match equi radii
    elseif strcmp(Pext_type,'ga')
        Pext_Amp_Freq = [P_guess; dt; tw;];
    end

    disptime = 0; % 1 = Displays time to complete simulation
    Dim = 0;  % 1 = displays results in dimensional form
    comp = 1; % 0 uses Rayleigh-Plesset, 1 uses Keller-Miksis

    if strcmp(Pext_type,'ga')
        Rmax = R0;
        R0 = eqR;
    end

    % Now, feed info to solver:
    %Variables:
    %t2 = simulation time
    %R2 = simulation radius
    %R2dot = velocity of bubble wall
    %P = bubble pressure
    %S = stress integral
    %T = temperature inside bubble
    %C = relative concentration of vapor
    %Tm = temperature of wall
    [t2, R2, R2dot, P, S, T, C, Tm, tdel, Tdel, Cdel] = IMRsolver(model, G, G1, mu,...
        tspan, R0, NT, NTM,  Pext_type,Pext_Amp_Freq , disptime, ...
        Tgrad, Tmgrad, Cgrad, Dim, comp, IX, RMesh); % ZZ - Add initial expansion, spatial mesh


expts(i).t_sims = t2;
expts(i).R_sims = R2 ;
expts(i).G_sims = G;
expts(i).mu_sims = mu ;
% Generate plots:

% Colors
pb = [43,131,186]/255; % Pretty blue
pr = [215,25,28]/255; % Pretty red
pg = [0,136,55]/255; % Pretty green
pp = [117,107,177]/255; % Pretty purple
lb = [107,174,214]/255;
lr = [251,106,74]/255;
lg = [153,216,201]/255;
lp = [188,189,220]/255;
%
% figure('units','normalized','outerposition',[0 0 0.6 0.7]);
% grid on; box on;
% 
% plot(t*1e6,R*1e6,'o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',pb,'Color',pb);hold on;
% plot(t2*1e6,R2*1e6,'-','LineWidth',2.5,'Color',pr) % Fit for B scale
% xlim([0,100]);
% ylim([0,350]);
% xl1 =  xlabel('Time $(\mu s)$','interpreter','latex');         
% yl1 = ylabel(' Raduis $(\mu m)$','interpreter','latex');
% set(gca,'FontSize',20);
% set(gca,'TickLabelInterpreter','latex');
% set(xl1,'Interpreter','Latex','FontSize',20,'FontWeight','bold');
% set(yl1,'Interpreter','Latex','FontSize',20,'FontWeight','bold');
% drawnow;
% lg0 = legend('Experiment','Simulation');
% set(lg0,'FontSize',20,'Location','northeast','Interpreter','Latex','Orientation','vertical');
% 
% disp([fileName,' is Done! Check results.']);

% savefig(fullfile(saveFigs,[fileName,'.fig']))
% saveas(gcf,fullfile(savePng,[fileName,'.png']))
end
% close all
%
hold on
% plot(t2*1e6/(23.05),R2/2.4475e-4)
% plot(t2/(R0/Uc),R2/R0,'rs')
plot(t2,R2,'rs')

%

% t_collapse_Plots;