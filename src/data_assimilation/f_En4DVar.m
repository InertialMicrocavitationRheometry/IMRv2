% file f_En4DVar.m
% brief file contains the Ensemble 4D variation Kalman filter

% brief This function conducts the IMR data assimilation for En4D Kalman
% filter code
function [x,ensemble,p_post,params]= f_En4DVar(param_pre,t,yth,R0_all, ...
    Req_all,tspan_all,peak_time_idx)

% number of ensemble members
q               =    48;
exp_i           =    1;

% data assimilation parameters

% data assimilation method ('En4D','EnKS',('EnKF'))
% method          =   'En4D';

% initial parameter guesses (all must be specified even if not used, in
% order to run):
model           =    param_pre{1};
P_prior         =    param_pre{3};
mu_theta        =    P_prior.mu;
sigma_theta     =    P_prior.sigma;
theta_params    =    mvnrnd(mu_theta,sigma_theta,q);
G_prior         =    theta_params(:,1) ;
mu_prior        =    theta_params(:,2)  ;
alpha_prior     =    theta_params(:,3);
G_guess         =    mean(G_prior);
mu_guess        =    mean(mu_prior);
alpha_guess     =    mean(alpha_prior);
Input_prior     =    true;
G1_guess        =    1e9;
lambda_nu_guess =    0.1;
% leave as 2, initializes ensemble with truth + noise
% init_scheme     =    2;

% ending criteria for iterative optimization:

% threshold difference in norm of state vector between steps
epsilon         =    1e-10;
% max # of iterations until end optimization (5 to 10 is usually good)
max_iter        =    5;

% Note: Beta coefficients are fixed as equal here. They can be modified in
% the main_En4D_peaks and main_mda_peaks files if needed, if more weight
% needs to be attributed to earlier or later points in assimilation

% modeling parameters

% amount of nodes inside the bubble
NT = 240;
% amount of nodes outside the bubble
NTM = 240;
% type of external forcing
Pext_type = 'IC';
% (N/m) Liquid Surface Tension
ST = 0.056;
% thermal effects inside bubble
Tgrad = 1;
% thermal effects outside bubble
Tmgrad = 1;
% vapor diffusion effects
Cgrad = 1;
% activates the effect of compressibility (0=Rayleigh-Plesset, 1=Keller-Miksis)
comp = 1;
% displays timesteps in En4D run (for debugging)
% disp_timesteps = 1;

% following should not be changed (untested):

% 1 = display simulation time
disptime = 0;
% 1 = output variables in dimensional form
% Dim = 0;

% covariance Inflation parameters
% the following is only used for the IEnKS. Default values provided below
% see Spratt et al. (2020) section 3.3 for details on theta, lambda

% 1 is scalar alpha CI, 2 is RTPS (BEST)
% CI_scheme = 2;
% multiplicative covariance parameter (0.5 < theta < 0.95)
CI_theta = 0.7;
% set to 1 for additive covariance (else 0)
% CI_add = 0;

% additive covariance parameter (lambda in paper) (1.005 < beta < 1.05)
beta = 1.02;

% random noise in forecast step (only for EnKF)
alpha = 0.005;

% spread of parameters in the ensemble
Rspread = 0.0;
Uspread = 0.0;
Pspread = 0.0;
Sspread = 0.0;
tauspread = 0.0;
Cspread = 0.0;
Tmspread = 0.0;
Brspread = 0.0;
Fohspread = 0.0;
% set to 0 if not used in model
Despread = 0;
% set to 0 if not used in model
lambda_nuspread = 0;

% do not modify
visco_params = struct('G',G_guess,'G1',G1_guess,'mu',mu_guess, ...
    'alpha',alpha_guess,'lambda_nu',lambda_nu_guess);
est_params = [];

% initialize and import data
clear aux1 aux2 dy2;
n =  length(t)-1;

% pre-allocating memory
x          =    zeros(2*NT+NTM+11,q,n+1);
x_est      =    zeros(2*NT+NTM+11,n+1);
E_est      =    zeros(2*NT+NTM+11,q,max_iter);
ensemble   =    zeros(q,max(peak_time_idx)-1);

tspan      =    tspan_all(exp_i);
R0         =    R0_all(exp_i);
Req        =    Req_all(exp_i);











% this file sets up an initial ensemble for experimental data setups

% state vector used: [R,U,P,S,Tau,C,Tm, 1/Ca, 1/Re]

addpath ./IMR-vanilla/functions
addpath ./IMR-vanilla/src

% shuffle rng to ensure randomness of results
rng('shuffle');

% guess for parameters

%NT = 500; %400 % amount of nodes inside the bubble (>=500 is a good to start)
%NTM = 500; %20 % amount of nodes outside the bubble (>=10 is a good to start)
%Pext_type = 'IC'; % Type of external forcing. Refer to RP_Cav

% find Req and calculate initial partial pressure needed parameters
% ST = 0.056; % (N/m) Liquid Surface Tension
Pmt_temp = IMRcall_parameters(R0,G_guess,G1_guess,mu_guess); % Calls parameters script
Ca = Pmt_temp(5);
Re = Pmt_temp(6);
P_inf = Pmt_temp(19);
T_inf = Pmt_temp(20);

% Req = mean(yth(end-5:end)); % take mean of end of sim to be Req
P_guess = (P_inf + (2*ST)/(Req*R0) - Pvsat(T_inf))*(Req^3);

Pext_Amp_Freq =[P_guess 0];
% [ Pressure ; Freq ], Freq = 0 for IC
% Pext_Amp_Freq = [100 0];













% determine initial state vector based on parameters
% This script determines initial X0 and x0
%{
global tspan R0 NT NTM Pext_type Pext_Amp_Freq Tgrad Cgrad model G G1 ...
    mu t0 neoHook nhzen sls linkv k chi Fom Foh We Br A_star B_star ...
    Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star De deltaY ...
    yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm x0_true N
%}
%***************************************
% Extract viscoelastic parameters from struct
G = visco_params.G;
G1 = visco_params.G1;
mu = visco_params.mu;
alpha = visco_params.alpha;
lambda_nu = visco_params.lambda_nu;

% Load Parameters :
Pmt = IMRcall_parameters(R0,G,G1,mu); % Calls parameters script
k = Pmt(1);
chi = Pmt(2);
Fom = Pmt(3);
Foh = Pmt(4);
Ca = Pmt(5);
Re = Pmt(6);
We = Pmt(7);
Br = Pmt(8);
A_star = Pmt(9);
B_star = Pmt(10);
Rv_star = Pmt(11);
Ra_star = Pmt(12);
P0_star = Pmt(13);
t0 = Pmt(14);
C0 = Pmt(15);
L = Pmt(16);
L_heat_star = Pmt(17);
Km_star = Pmt(18);
P_inf = Pmt(19);
T_inf = Pmt(20);
C_star = Pmt(21);
De = Pmt(22);
rho = Pmt(23);

%****************************************

% Material Choice
neoHook = 0;
nhzen = 0;
sls = 0;
linkv = 0;
fung = 0;
fung2 = 0;
fungexp = 0;
fungnlvis = 0;
if strcmp(model,'neoHook') == 1
    neoHook = 1;
elseif strcmp(model,'fung') == 1
    fung = 1;
elseif strcmp(model,'fung2') == 1
    fung2 = 1;
elseif strcmp(model,'fungexp') == 1
    fungexp = 1;
elseif strcmp(model,'fungnlvis') == 1
    fungnlvis = 1;
elseif strcmp(model,'nhzen') == 1
    nhzen = 1;
elseif strcmp(model,'sls') == 1
    sls = 1;
elseif strcmp(model,'linkv') == 1
    linkv = 1;
else
    nhzen = 1;
end


% Needed to account for fast diffusion
P0_star = P0_star - (1-Cgrad)*Pvsat(1*T_inf)/P_inf;

% When we assume water vapor undergoes infinitely fast mass diffusion
% the vapor pressure is constant and P is the pressure of
% non-condesible gas

%******************************************
% Creates finite difference matrices
D_Matrix_T_C = Finite_diff_mat(NT,1,0);
DD_Matrix_T_C = Finite_diff_mat(NT,2,0);
D_Matrix_Tm = Finite_diff_mat(NTM,1,1);
DD_Matrix_Tm = Finite_diff_mat(NTM,2,1);
%******************************************

%******************************************
% Create spatial nodes

% Inside the bubble
N = NT-1;
deltaY = 1/N;
i = 1:1:N+1;
yk = ((i-1)*deltaY)';

% Outside the bubble
Nm = NTM-1;
deltaYm = -2/Nm;
j = 1:1:Nm+1;
xk = (1+(j-1)*deltaYm)';
yk2 = ((2./(xk+1)-1)*L+1);

%******************************************
% Initial Conditions
tspan_star = tspan/t0;
R0_star = 1;
U0_star = 0;  % Change as needed
%Z10 = 0;
S0 = 0;
Tau0 = zeros(1,NT);
C0 = C0*ones(1,NT);
Tm0 = ones(1,NTM);
if strcmp(Pext_type,'ga')
    dt_star = Pext_Amp_Freq(2)/t0;
    w_star = Pext_Amp_Freq(3)*t0;
end

% Need to modify initial conditions for the Out-of-Equilibrium Rayleigh
% Collapse:
if strcmp(Pext_type,'IC')
    Pv = Pvsat(1*T_inf)/P_inf;
    P0_star = Pext_Amp_Freq(1)/P_inf + Cgrad*Pvsat(1*T_inf)/P_inf;
    % Need to recalculate initial concentration
    theta = Rv_star/Ra_star*(P0_star-Pv)/Pv; % mass air / mass vapor
    C0 = 1/(1+theta);
    
    % Calculate the equilibrium radii ratio for initial stress state:
    if Req == 0
        [REq,~,~] = IMRCalc_Req(R0, Tgrad, Cgrad, Pext_Amp_Freq(1), G, G1, mu);
    end
    REq = Req;
    %REq = 1; %removed 6/15/16 by Jon
    C0 = C0*ones(1,NT);
    %U0_star = -1*(1-P0_star)/(C_star); %Initial velocity due to shockwave
    U0_star = 0;
    
    if sls == 1 || linkv == 1
        S0 = -4/(3*Ca)*(1-REq^3);
    elseif nhzen == 1 || neoHook == 1
        S0 = -1/(2*Ca)*(5-REq^4-4*REq);
    elseif fung == 1 || fungnlvis == 1
        S0 = -(1-3*alpha)*(5 - 4*REq - REq^4)/(2*Ca) - ...
            2*alpha*(-27/40 - 1/8*REq^8 - 1/5*REq^5 -1*REq^2 + 2/REq)/(Ca);
    elseif fung2 == 1
        S0 = 2*(1-3*alpha+4.5*alpha^2)*(-5/4 + 1*REq + 1/4*REq^4)/(Ca) + ...
            2*(27/40*alpha-221/90*alpha^2 + alpha^2/24*REq^12 + alpha^2/18*REq^9 + (alpha-3*alpha^2)/8*REq^8 + ... )
        2*alpha^2/6*REq^6 + (alpha-3*alpha^2)/5*REq^5 + 2*alpha^2/3*REq^3 + (2*alpha-6*alpha^2)/2*REq^2 - ...
            2*alpha^2*log(1/REq) - (2*alpha-6*alpha^2)/REq - 2/3*alpha^2/REq^3)/(Ca);
    elseif fungexp == 1
        tempbeta = linspace(1/REq,1,1e5);
        tempS = 2*(tempbeta.^-5 + tempbeta.^-2) .* exp(alpha*(tempbeta.^-4+2*tempbeta.^2-3))/(Ca);
        % figure, plot(tempbeta,tempS);
        S0 = (trapz(tempbeta,tempS));
    end
    
end

X0 = [R0_star U0_star P0_star S0 Tau0 C0 Tm0];

tau_del = [];

x0_true = [X0,Br,Foh,G,mu,De,alpha,lambda_nu,est_params];
N = length(x0_true);














% create initial ensemble

% custom spread:
x_init = x0_true';

if Input_prior
    Caspread    = 0;
    Respread    = 0;
    alphaspread = 0;
end

spread = [Rspread;
Uspread;
Pspread;
Sspread;
ones(NT,1)*tauspread; ...
    ones(NT,1)*Cspread;
ones(NTM,1)*Tmspread;
Brspread;
Fohspread; ...
    Caspread;
Respread;
Despread;
alphaspread;
lambda_nuspread];

xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
    repmat([0;
0;
0;
0;
zeros(2*NT+NTM,1);
0;
0;
0;
0;
0;
0;0],1,q) .* randn(N,q);

% using truncated distribution
if Input_prior
    xi(2*NT+NTM+7,:)      =   G_prior';
    xi(2*NT+NTM+8,:)      =   mu_prior';
    xi(2*NT+NTM+10,:)     =   alpha_prior';
else
    pd_G                  =   makedist('Normal','mu',G_guess,'sigma',G_guess*Caspread);
    t_G                   =   truncate(pd_G,1e-10,inf);                                  % truncated Gaussian distribution
    
    if Caspread ==0
        xi(2*NT+NTM+7,:)  =   G_guess;
    else
        xi(2*NT+NTM+7,:)  =   random(t_G,[1, q]);
    end
    
    pd_mu                 =   makedist('Normal','mu',mu_guess,'sigma',mu_guess*Respread);
    t_mu                  =   truncate(pd_mu,1e-10,inf);
    xi(2*NT+NTM+8,:)      =   random(t_mu,[1, q]);
    pd_alpha              =   makedist('Normal','mu',alpha_guess,'sigma',alpha_guess*alphaspread);
    
    if alphaspread ==0
        xi(2*NT+NTM+10,:) = alpha_guess;
    else
        xi(2*NT+NTM+10,:) = random(t_alpha,[1, q]);
    end
    
end

xi(3,:) = log(xi(3,:));

% constrain the pressure to be positive
xi(3,:) = log(x0_true(3));

Uc      = sqrt(P_inf/rho);

xi      = [xi(1:2*NT+NTM+6,:);...
    (xi(2*NT+NTM+7,:))./P_inf; ...
    ((xi(2*NT+NTM+8,:)).*Uc)./(P_inf*R0);...
    xi(end-2,:);...
    (xi(end-1,:));...
    xi(end,:)]; % for now












tspan_star =    max(tspan);
x(:,:,1)   =    xi;
tau_del    =    cell(q,1);
idx        =    1;

% iterate
vars = {NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
    Cgrad comp t0 neoHook nhzen sls linkv k chi Fom Foh We Br A_star ...
    B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
    De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm tspan_star NTM rho R0 fung fung2 fungexp fungnlvis};

opts.POSDEF = true;
opts.SYM = true;


%delta = 1.005; %inflation param
delta = beta;
%epsilon = 0.001;

% size of data assimilation window
l                   =   max(peak_time_idx(:)) - 1;
% MDA coefs (equal here)
Beta                =   (1/l).*ones(1,l);
timestep_time       =   zeros(1,n-l+1);
time_index          =   1;
% numbers of DA cycles
j_max               =   3;
params              =   [];

for j = 1:j_max
    
    % restrict mu and G to be positive
    if min(x(2*NT+NTM+7,:,j))<=0
        Idx = (x(2*NT+NTM+7,:,j))<=0;
        x(2*NT+NTM+7,Idx,j) = x(2*NT+NTM+7,Idx,j)-min(x(2*NT+NTM+7,:,j))+1e-6;
    end
    if min(x(2*NT+NTM+8,:,j))<=0
        Idx = (x(2*NT+NTM+8,:,j))<=0;
        x(2*NT+NTM+8,Idx,j) = x(2*NT+NTM+8,Idx,j)-min(x(2*NT+NTM+8,:,j))+1e-6;
    end
    
    %j print index for debugging if it crashes
    %timestep_time(j) = toc;
    x10 =  mean(squeeze(x(:,:,j)),2);
    A10 =  squeeze(x(:,:,j)) - x10*ones(1,q);
    x1  =  x10;
    TT  =  eye(q);
    dx1 =  1; %initialize greater than epsilon
    jj  =  1;
    A1  =  A10*TT;
    E1  =  x1*ones(1,q) + A1;
    
    for j1 = 1:size(A1,1)
        if std(A1(j1,:)) ==0
            std_prior(j1) =0;
        else
            std_prior(j1) = sqrt(robustcov(A1(j1,:)));
        end
    end
    
    while norm(dx1) > epsilon && jj <= max_iter
        % saving parameters at each iteration
        G_all            =  E1(2*NT+NTM+7,:)*P_inf;
        mu_all           =  E1(2*NT+NTM+8,:)/Uc*(P_inf*R0);
        alpha_all        =  (E1(2*NT+NTM+10,:));
        params(:,1,idx)  =  G_all;
        params(:,2,idx)  =  mu_all;
        params(:,3,idx)  =  alpha_all;
        idx              =  idx+1;
        
        TTinv = linsolve(TT,eye(size(TT)));
        t1    = t(exp_i,time_index);
        t2    = t(exp_i,time_index+l);
        
        parfor memb = 1:q
            [t_memb{memb}, EE{memb},~] =  f_new(t1,t2,E1(:,memb),vars,tau_del{memb});
            t_sim                      =  t_memb{memb};
            y_sim                      =  EE{memb}(:,1);
            U_sim                      =  EE{memb}(:,2);
            y2(1,memb,:)               =  interp1(t_sim,y_sim, t(exp_i,(1:l)+time_index), 'makima' );
        end
        
        clear y2b dy2 HA2
        
        for kk = 1:l
            y2b(:,kk)      =  mean(y2(:,:,kk),2);
            HA2(:,:,kk)    =  y2(:,:,kk) - y2b(:,kk)*ones(1,q);
            HA2(:,:,kk)    =  HA2(:,:,kk)*TTinv;
            dy2(:,kk)      =  mean(yth(:,time_index+kk)) - y2b(:,kk);
        end
        
        ensemble(:,1:l)           =  squeeze(y2);
        aux2                =  zeros(q,l);
        
        for kk = 1:l
            [R(kk),~]      =  robustcov(yth(:,time_index+kk));
            Rinv           =  linsolve(R(kk),eye(size(R(kk))));
            aux1(:,:,kk)   =  (HA2(:,:,kk)'*Beta(kk)*Rinv*HA2(:,:,kk))/(q-1);
            aux2(:,kk)     =  aux2(:,kk) + (HA2(:,:,kk)'*Beta(kk)*Rinv*(dy2(:,kk)))/(q-1);
        end
        
        GGinv               =  eye(q) + sum(aux1,3);
        GG                  =  linsolve(GGinv,eye(size(GGinv)));
        b                   =  linsolve(GGinv,sum(aux2,2));
        dx1                 =  A10*b + A10*GG*pinv(A10'*A10)*A10'*(x10-x1);
        x1                  =  x1 + dx1;
        TT                  =  sqrtm(GG);
        
        x_est(:,jj)         =  x1;
        timestep_time(jj)   =  toc;
        E_est(:,:,jj)       =  E1;
        A1                  =  A10*TT;
        E1                  =  x1*ones(1,q) + A1;
        
        if min(E1(2*NT+NTM+8,:))<=0
            Idx = (E1(2*NT+NTM+8,:))<=0;
            E1(2*NT+NTM+8,Idx) = E1(2*NT+NTM+8,Idx)-min(E1(2*NT+NTM+8,:))+1e-6;
        end
        disp(['ended iteration ',num2str(jj),' with norm(dx1) = ', ...
            num2str(norm(dx1)), ' and norm(dy2) = ' num2str(norm(dy2)), ' at ',num2str(timestep_time(jj)),' seconds for the ' num2str(exp_i) 'th experiments '])
        jj = jj+1;
    end
    
    % perform covaraince inflation:
    if j<j_max
        for j1 = 1:size(A1,1)
            if std(A1(j1,:)) ==0
                std_prior(j1) =0;
            else
                std_post(j1) = sqrt(robustcov(A1(j1,:)));
            end
        end
        RTPS              = 1+CI_theta*(std_prior-std_post)./std_post;
        RTPS(isnan(RTPS)) = 1;
        RTPS(isinf(RTPS)) = 1;
        for j1 = 1:size(A1,1)
            A1(j1,:) = mvnrnd(0,(RTPS(j1).*std_post(j1))^2,q);
        end
        Input_prior =  false;
        E1 = x1*ones(1,q) + A1;
    end
    x(:,:,j+1) = E1;
end

% posterior ensembles:
G_post = E1(2*NT+NTM+7,:)*P_inf;
mu_post = E1(2*NT+NTM+8,:)/Uc*(P_inf*R0);
alpha_post = (E1(2*NT+NTM+10,:));
[sigma_post,mu_post] = robustcov([G_post;
mu_post;
alpha_post]');
p_post.mu = mu_post;
p_post.sigma = sigma_post;

end
