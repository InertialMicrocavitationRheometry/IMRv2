% equation options
radial          = 1;        % 1 : RP, 2 : K-M P, 3 : K-M E, 4 : Gil
bubtherm        = 0;        % 0 : polytropic assumption, 1: thermal PDE in bubble
medtherm        = 0;        % 0 : cold fluid, 1: warm fluid assumption
stress          = 0;        % 1 : NHKV, qKV, 2: linear Maxwell, Jeffreys, Zener, 3: UCM or OldB, 4: PTT, 5: Giesekus
eps3            = 0;        % this value must be (0, 0.5]
vapor           = 0;        % 0 : ignore vapor pressure, 1 : vapor pressure
masstrans       = 0;        % mass transfer, default is no mass transfer 
    
% solver options
TFin            = 20e-6;     % final time (s)
TVector         = [0 TFin];
method          = 45;       % ode45 setting for the time stepper
spectral        = 0;        % force spectral collocation solution
divisions       = 0;        % minimum number of timesteps
Nt              = 10;       % number of points in bubble, thermal PDE
Mt              = 10;       % number of points outside of bubble, thermal PDE
Nv              = 150;      % number of points outside of bubble, viscoelastic stress PDE     
Lv              = 3;        % characteristic length for grid stretching, leave at 3
Lt              = 3;        % characteristic length for grid stretching, leave at 3
    
% initial conditions
R0              = 100E-6;   % initial bubble radius
U0              = 0;        % initial velocity (m/s)
Req             = 50E-6;    % Equilibrium radius for pre-stress bubble, see Estrada JMPS 2017    

% output options
dimensionalout  = 0;        % output result in dimensional variables
progdisplay     = 0;        % display progress while code running

% acoustic options
rho8            = 1064;%997;              % far-field density (kg/m^3)   
GAM             = 3049.13*1e5;      % state equation parameter (Pa)
nstate          = 7.15;             % state equation parameter
P8              = 101325;           % far-field pressure (Pa)
C8              = 1484;             %sqrt(nstate*(P8 + GAM)/rho8); % far-field sound speed (m/s)

% pressure wave options
pA              = 0*1e6;      % pressure amplitude (Pa)
omega           = 0*4e6*2*pi; % frequency (rad/s)
TW              = 0;        % gaussian width (s)
DT              = 0;        % delay (s)
mn              = 0;        % power shift for waveform
wave_type       = 2;        % wave type oscillating bubble, see f_pinfinity
    
% stress options
S               = 0.072;%0.056;% 0.072              % (N/m) Liquid Surface Tension 
vmaterial       = 'water';
G               = 312.5;%1E3;                % (Pa) Medium Shear Modulus 
lambda1         = 0*0.5e-5;             % relaxation time (s)
lambda2         = 0;            % retardation time (s)
mu8             = 0.027606;%0.0246;%1E-3;
alphax          = 0;        % qKV term
S0              = 0;
%(3*alphax-1)*(5 - (Req/1)^4 - 4*(Req/1))/(2*Ca) + ...
% (2*alphax/Ca)*(27/40 + (1/8)*(Req/1)^8 + (1/5)*(Req/1)^5 + (1/2)*(Req/1)^2 - ...
% 2*1/Req)
    
% thermal options

% gas properties
kappa           = 1.4;               % Ratio of Specific Heats 
% medium properties 
AT              = 5.28e-5;            % (W/m-K^2)Thermal Conductivity coeff
BT              = 1.165e-2;           % (W/m-K)Thermal Conductivity coeff
T8              = 298.15;             % (K) Far field temp. 
Km              = 0.55;%0.615;        % (W/m-K)Thermal Conductivity Medium
Cp              = 4181;               % Specific Heat Medium J/Kg K;
Dm              = Km / (rho8*Cp) ;    % Thermal Diffusivity m^2/s 
	
% mass transfer options
D0              = 24.2e-6;            % Diffusion Coeff m^2/s
L_heat          = 2264.76e3;          % (J/Kg) Latent heat of evaporation
Ru              = 8.3144598;          % (J/mol-K) Universal Gas Constant % Ru/(18.01528e-3);%Ru/(28.966e-3);     
Rv              = Ru/(18.01528e-3);   % (J/Kg-K) Gas constant vapor
Ra              = 438.275;            % (J/Kg-K)Gas constant air    