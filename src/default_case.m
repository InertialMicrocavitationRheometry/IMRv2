% equation options

% 1 : RP, 2 : K-M P, 3 : K-M E, 4 : Gil
radial          = 1;
% 0 : polytropic assumption, 1: thermal PDE in bubble
bubtherm        = 0;
% 0 : cold fluid, 1: warm fluid assumption
medtherm        = 0;
% 1 : NHKV, qKV, 2: linear Maxwell, Jeffreys, Zener, 3: UCM or OldB, 4: PTT, 5: Giesekus
stress          = 0;
% this value must be (0, 0.5]
eps3            = 0;
% 0 : ignore vapor pressure, 1 : vapor pressure
vapor           = 0;
% mass transfer, default is no mass transfer
masstrans       = 0;

% solver options

% final time (s)
TFin            = 20e-6;
% time span vector
TVector         = [0 TFin];
% ode45 setting for the time stepper
method          = 23;
% force spectral collocation solution
spectral        = 0;
% minimum number of timesteps
divisions       = 0;
% number of points in bubble, thermal PDE
Nt              = 25;
% number of points outside of bubble, thermal PDE
Mt              = 25;
% number of points outside of bubble, viscoelastic stress PDE
Nv              = 150;
% characteristic length for grid stretching, leave at 3
Lv              = 3;
% characteristic length for grid stretching, leave at 3
Lt              = 2;

% initial conditions

% initial bubble radius
R0              = 100E-6;
% initial velocity (m/s)
U0              = 0;
% Equilibrium radius for pre-stress bubble, see Estrada JMPS 2017
Req             = 50E-6;

% output options

% output result in dimensional variables
dimensionalout  = 0;
% display progress while code running
progdisplay     = 0;

% acoustic options

% far-field density (kg/m^3)
rho8            = 1000;
% state equation parameter (Pa)
GAM             = 3049.13*1e5;
% state equation parameter
nstate          = 7.15;
% far-field pressure (Pa)
P8              = 101325;
% sqrt(nstate*(P8 + GAM)/rho8); % far-field sound speed (m/s)
C8              = 1484;

% pressure wave options

% pressure amplitude (Pa)
pA              = 0*1e6;
% frequency (rad/s)
omega           = 0*4e6*2*pi;
% Gaussian width (s)
TW              = 0;
% delay (s)
DT              = 0;
% power shift for waveform
mn              = 0;
% wave type oscillating bubble, see f_pinfinity
wave_type       = 2;

% stress options

% (N/m) Liquid Surface Tension
S               = 0.072;
vmaterial       = 'water';
% (Pa) Medium Shear Modulus
G               = 1e3;
% relaxation time (s)
lambda1         = 1e-8;
% retardation time (s)
lambda2         = 1e-8;
% qKV term
alphax          = 0.25;

% thermal options

% gas properties

% ratio of specific heats
kappa           = 1.47;

% medium properties

% (W/m-K^2) thermal conductivity coeff
AT              = 5.28e-5;
% (W/m-K) thermal conductivity coeff
BT              = 1.165e-2;
% (K) far field temperature
T8              = 300;
% (W/m-K) thermal conductivity medium
Km              = 0.615;
% (J/Kg K) specific heat medium
Cp              = 3.61e3;
% (m^2/s) thermal diffusivity
Dm              = Km / (rho8*Cp);

% mass transfer options

% (m^2/s) diffusion coefficient
D0              = 24.2e-6;
% (J/Kg) latent heat of evaporation
L_heat          = 2264.76e3;
% (J/mol-K) universal gas constant
Ru              = 8.3144598;
% (J/Kg-K) gas constant vapor
Rv              = Ru/(18.01528e-3);
% (J/Kg-K) gas constant air
Ra              = 438.275;

% viscosity variables
[mu8,Dmu,v_a,v_nc,v_lambda,vmat] = f_nonNewtonian_Re(vmaterial);

% pressure variables
Pv              = f_pvsat(T8);
P0              = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
