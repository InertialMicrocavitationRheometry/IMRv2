% file f_call_params.m
% brief contains function f_call_params

% brief This function sets up the simulation environment. The function can
% receive an input file. The input file must be in the formatted similar to
% the default_case.m file. Otherwise, if no input file is provided, the
% code will read in the default case. The default case can be altered based
% on the input arguments when invoking either the finite difference or
% spectral module.
function [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
    acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] = ...
    f_call_params(varargin)

disp('--- Inertial Microcavitation Rheometry forward solver ---');
% check that all inputs are matched
if mod(nargin,2) == 1
    error('input error: unmatched inputs');
end

% checking for the casefile, if any
defaultread = true;
for n = 1:2:nargin
    if strcmpi(varargin{n},'casefile') == 1
        cfname = varargin{n+1};
        disp('Case file: Using given casefile',cfname);
        try
            run(cfname);
        catch
            error('failed to read the given file');
        end
        defaultread = false;
        break;
    end
end

% otherwise, reading the default casefile
if defaultread
    disp('Case file: Using default case file');
    run('default_case.m');
end

% overrides defaults with options and dimensional inputs %

% load inputs
misscount = 0;
p0set = 0;
tempset = 0;
dmset = 0;
tflag = 0;
visflag = 0;
for n = 1:2:nargin
    if strcmpi(varargin{n},'p0') == 1
        p0set = 0;
    end
    switch lower(varargin{n})
        % equation options
        case 'radial',      radial = varargin{n+1};
        case 'bubtherm',    bubtherm = varargin{n+1};
        case 'medtherm',    medtherm = varargin{n+1};
        case 'stress',      stress = varargin{n+1};
        case 'eps3',        eps3 = varargin{n+1};
        case 'vapor',       vapor = varargin{n+1};
        case 'masstrans',   masstrans = varargin{n+1};
        
        % solver options
        case 'method',      method = varargin{n+1};
        case 'spectral',    spectral = varargin{n+1};
        case 'divisions',   divisions = varargin{n+1};
        case 'nv',          Nv = varargin{n+1};
        case 'nt',          Nt = varargin{n+1};
        case 'mt',          Mt = varargin{n+1};
        case 'lv',          Lv = varargin{n+1};
        case 'lt',          Lt = varargin{n+1};
        case 'tfin',        TFin = varargin{n+1};
        tflag = tflag + 1;
        case 'tvector',     TVector = varargin{n+1};
        tflag = tflag + 1;
        TFin = 0;
        
        % initial options
        case 'collapse',    collapse = varargin{n+1};
        case 'r0',          R0 = varargin{n+1};
        P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
        case 'u0',          U0 = varargin{n+1};
        case 'req',         Req = varargin{n+1};
        P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
        case 'stress0',     Szero = varargin{n+1};
        
        % output options
        case 'dimout',      dimensionalout = varargin{n+1};
        case 'progdisplay', progdisplay = varargin{n+1};
        
        % acoustic options
        case 'rho8',        rho8 = varargin{n+1};
        case 'gam',         GAM = varargin{n+1};
        case 'nstate',      nstate = varargin{n+1};
        case 'p8',          P8 = varargin{n+1};
        case 'c8',          C8 = varargin{n+1};
        
        % pressure waveform options
        case 'pa',          pA = varargin{n+1};
        case 'omega',       omega = varargin{n+1};
        case 'tw',          TW = varargin{n+1};
        case 'dt',          DT = varargin{n+1};
        case 'mn',          mn = varargin{n+1};
        case 'wave_type',   wave_type = varargin{n+1};
        
        % stress options
        case 'mu',          mu8 = varargin{n+1};
        visflag = visflag + 1;
        case 'g',           G = varargin{n+1};
        case 'lambda1',     lambda1 = varargin{n+1};
        case 'lambda2',     lambda2 = varargin{n+1};
        case 'alphax',      alphax = varargin{n+1};
        case 'surft',       S = varargin{n+1};
        % P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
        
        % viscosity options
        case 'du',          Dmu         = varargin{n+1};
        case 'mu0',         muo         = varargin{n+1};
        Dmu = muo - mu8;
        case 'v_a',         v_a         = varargin{n+1};
        case 'v_nc',        v_nc        = varargin{n+1};
        case 'v_lambda',    v_lambda    = varargin{n+1};
        case 'nu_model',    nu_model     = varargin{n+1};
        
        % thermal options
        case 't8',          T8 = varargin{n+1};
        tempset = 1;
        case 'kappa',       kappa = varargin{n+1};
        case 'atg',         ATg = varargin{n+1};
        case 'btg',         BTg = varargin{n+1};
        case 'atv',         ATv = varargin{n+1};
        case 'btv',         BTv = varargin{n+1};
        case 'km',          Km = varargin{n+1};
        case 'dm',          Dm = varargin{n+1};
        dmset = 1;
        
        % mass transfer options
        case 'dmass',       D0 = varargin{n+1};
        case 'lheat',       L_heat = varargin{n+1};
        case 'rv',          Rv = varargin{n+1};
        case 'ra',          Ra = varargin{n+1};
        
        % pressure options
        case 'pv',          Pv = varargin{n+1};
        P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
        case 'p0',          P0 = varargin{n+1};
        
        otherwise,          misscount = misscount + 1;
        
    end
end

if tempset == 1
    % recalculating the vapor pressure
    Pv = vapor*f_pvsat(T8);
    P0 = (P8 + 2*S/Req - Pv)*(Req/R0)^3;
end
if p0set == 0
    % need to add Pv_sat at room temp
    Pv = vapor*f_pvsat(T8);
    P0 = (P8 + 2*S/Req - Pv)*(Req/R0)^3;
end
if dmset == 0
    % (m^2/s) thermal diffusivity
    Dm = Km / (rho8*Cp);
end

% loading waveform data
if wave_type < 0
    if wave_type == -1
        waveform_dir = './waveform_data/hn.mat';
    elseif wave_type == -2
        waveform_dir = './waveform_data/mn.mat';
    elseif wave_type == -3
        waveform_dir = './waveform_data/ml.mat';
    end
    wave_poly = load(waveform_dir,'poly');
    wave_dpoly = load(waveform_dir,'dpoly');
else
    wave_poly = [];
    wave_dpoly = [];
end

check = 1-isnumeric(radial);
if check || radial > 4 || radial <= 0
    error('INPUT ERROR: radial must be 1, 2, 3, or 4');
end
check = 1-isnumeric(bubtherm);
if check || bubtherm ~= 0 && bubtherm ~= 1
    error('INPUT ERROR: bubtherm must be 0 or 1');
end
check = 1-isnumeric(medtherm);
if check || medtherm ~= 0 && medtherm ~= 1
    error('INPUT ERROR: medtherm must be 0 or 1');
end
check = 1-isnumeric(stress);
if check || stress > 5 || stress < 0
    error('INPUT ERROR: stress must be between 0 to 5');
end
check = 1-isnumeric(vapor);
if check || vapor ~= 0 && vapor ~= 1
    error('INPUT ERROR: vapor must be 0 or 1');
end
check = 1-isnumeric(masstrans);
if check || masstrans ~= 0 && masstrans ~= 1
    error('INPUT ERROR: masstrans must be 0 or 1');
elseif check || masstrans == 1 && vapor == 0
    error('INPUT ERROR: if masstrans is 1 you must have vapor be 1')
end
check = 1-isnumeric(wave_type);
if check || wave_type > 5 || wave_type <= -3
    error('INPUT ERROR: wavetype must be between -3 and 5');
end
if (tflag > 1)
    error('INPUT ERROR: Only tvector or tfin can be specified, not both');
end
if TVector == 0
    TVector = [0 TFin];
end

% intermediate calculated variables

% far-field thermal conductivity (W/(m K))
K8      = ATg*T8+BTg;
% dimensional parameter for gas constants
Rnondim = P8/(rho8*T8);
% characteristic velocity (m/s)
Uc      = sqrt(P8/rho8);

% Final non-dimensional variables
Pref    = P8;
% dimensionless vapor and infinity pressure
Pv_star = vapor*Pv/Pref;
P0_star = P0/Pref;
% characteristic time (s)
t0      = R0/Uc;

% dimensionless waveform parameters
tvector = TVector./t0;
% non-dimensional frequency
om      = omega*t0;
ee      = pA/Pref;
tw      = TW*t0;
dt      = DT/t0;
% acoustic properties

% bulk liquid stiffness
GAMa    = GAM/P8;
% speed of sound
Cstar   = C8/Uc;
% thermal properties
chi     = T8*K8/(P8*R0*Uc);
iota    = Km/(K8*Lt);
Foh     = Dm/(Uc*R0);
alpha_g   = ATg*T8/K8;
beta_g    = BTg/K8;
alpha_v   = ATv*T8/K8;
beta_v   = BTv/K8;
Br      = Uc^2/(Cp*T8);
% mass diffusion
Fom     = D0/(Uc*R0);
% mass of vapor
mv0     = Pv*vapor*(4/3*pi*R0^3)/Rv/T8;
% mass of non-condensible gas
ma0     = P0*(4/3*pi*R0^3)/Ra/T8;
Mnondim = rho8*(4/3*pi*R0^3);
mv0     = mv0/Mnondim;
ma0     = ma0/Mnondim;
Rv_star = Rv/Rnondim;
Ra_star = Ra/Rnondim;
% mass air / mass vapor
theta = Rv_star/Ra_star*(P0_star-Pv_star)/Pv_star;
% initial vapor concentration
C0 = 1/(1+theta);
L_heat_star = L_heat/(Uc)^2;

% viscoelastic properties

% Cauchy number
Ca      = Pref/G;
% Reynolds number
Re8     = Pref*R0/(mu8*Uc);
if Dmu ~= 0
    DRe = Pref*R0/(Dmu*Uc);
else
    DRe = 0;
end
if DRe==0
    iDRe = 0;
else
    iDRe = 1/DRe;
end

% Weber number
log_We = log(0.5) + log(Pref) + log(R0) - log(S);
We = exp(log_We);
% relaxation time
v_lambda_star = v_lambda/t0;
% Weissenberg number
LAM     = lambda2/lambda1;
% Deborah number
De      = lambda1*Uc/R0;
% dimensionless initial conditions
Rzero   = 1;
Uzero   = U0/Uc;

% overwrite defaults with nondimensional inputs
if isempty(varargin) == 0
    for n = 1:2:nargin
        switch lower(varargin{n})
            % dimensionless state variables
            case 'cstar',   Cstar = varargin{n+1};
            case 'gama',    GAMa = varargin{n+1};
            % dimensionless waveform parameters
            case 'om',      om = varargin{n+1};
            case 'ee',      ee = varargin{n+1};
            case 'twx',     tw = varargin{n+1};
            case 'dtx',     dt = varargin{n+1};
            % dimensionless numbers
            case 'we',      We = varargin{n+1};
            case 're',      Re8 = varargin{n+1};
            case 'dre',     DRe = varargin{n+1};
            case 'ca',      Ca = varargin{n+1};
            case 'lam',     LAM = varargin{n+1};
            case 'de',      De = varargin{n+1};
            case 'foh',     Foh = varargin{n+1};
            case 'br',      Br = varargin{n+1};
            case 'fom',     Fom = varargin{n+1};
            % dimensionless thermal quantities
            case 'alpha_g', alpha_g = varargin{n+1};
            case 'beta_g',  beta_g = varargin{n+1};
            case 'alpha_v', alpha_v = varargin{n+1};
            case 'beta_v',  beta_v = varargin{n+1};
            case 'chi',     chi = varargin{n+1};
            case 'iota',    iota = varargin{n+1};
            % dimensionless initial conditions
            case 'rzero',   Rzero = varargin{n+1};
            case 'uzero',   Uzero = varargin{n+1};
            case 'p0star',  P0_star = varargin{n+1};
            otherwise
            misscount = misscount + 1;
        end
    end
end
% check that all inputs were accounted for
if misscount ~= nargin/2
    error(['INPUT ERROR: ' num2str(misscount-nargin/2) ' unrecognized input(s)']);
end

% final setting adjustments
if bubtherm == 0
    medtherm = 0;
end
if medtherm == 1
    bubtherm = 1;
end

% 1 : N-H, 2: qN-H, 3: linear Maxwell, Jeffreys, Zener, 5: UCM or OldB, 6: PTT, 7: Giesekus
if stress == 0 || stress == 1 || stress == 2 || stress == 3 || stress == 4
    spectral = 0;
    if stress == 3 || stress == 4
        Nv = 1;
    else
        Nv = 0;
    end
elseif stress == 5
    Nv = (Nv - 1)*(spectral == 1) + 1;
    Ca = -1;
elseif stress == 6
    Nv = (Nv - 1)*(spectral == 1) + 1;
    Ca = -1;
    spectral = 1;
elseif stress == 7
    Nv = (Nv - 1)*(spectral == 1) + 1;
    Ca = -1;
    spectral = 1;
end

if Ca == -1
    Ca = Inf;
end

if stress == 1 || stress == 2 || stress == 3 || stress == 4
    JdotA = 4/Re8;
elseif stress == 5 || stress == 6
    JdotA = 4*LAM/Re8;
else
    JdotA = 0;
end
if spectral == 1
    JdotA = 0;
end

% Keller-Miksis equation
% if linkv==1 || neoHook==1 || Yeoh==1
% JdotA = 4/Re8;
% elseif sls==1 || nhzen==1 || fdkv==1 || zzzen==1 || fdmax==1  % ZZ - Not exactly true for FDKV, but I also don't think this matters ...
    % JdotA = 0;
% elseif nhkv_pld==1
%     %JdotA = 4/Re8*(2^alpha_g+1)/3*(abs(Rdot)/R)^(alpha_g-1);
%     JdotA = 4/Re8/3*(2^alpha_g+1)*sign(Rdot)*(abs(Rdot)/R)^(alpha_g)*R^2/Rdot^2;
%     if isnan(SdotA)
%         SdotA=4/Re8;
%     end
% end

if stress == 0 || stress == 1 || stress == 2 || stress == 3 || stress == 4
    zeNO = 0;
else
    zeNO = 1;
end

if (stress == 3 || stress == 4) && De == 0
    error('INPUT ERROR: De can not equal zero for stress = 3 or 4');
end

% inertial Rayleigh collapse out of equilibrium initial conditions
opts = optimset('display','off');
if collapse && vapor
    
    Pv = f_pvsat(T8)/P8;
    % Calculate the equilibrium radii ratio for initial stress state:
    ma0 = P0_star/Ra_star;
    % parameterized function
    fun = @(x) Pv*(1+(ma0/x)*(Ra_star/Rv_star))-1-1/We*(Pv/(Rv_star*x))^(1/3);
    MTotal0 = P0_star/Ra_star + Pv/Rv_star;
    MVE = fzero(fun,MTotal0,opts);
    while (isnan(MVE))
        MTotal0 = MTotal0/1.11;
        MVE = fzero(fun,MTotal0,opts);
    end
    Pb_star = P0_star + Pv;
    % Need to recalculate initial concentration, mass air / mass vapor
    theta = Rv_star/Ra_star*(Pb_star/Pv-1);
    C0 = 1/(1+theta);
    Req = (Rv_star*MVE/Pv)^(1/3);
    
else
    
    Pv = f_pvsat(T8)/P8*vapor;
    Pb_star = P0_star + Pv;
    Req = Req/R0;
    
end

Pv_star = Pv;
Req_zero = Req;

% initial stress field for bubble collapse
if stress < 3
    Szero = [];
elseif stress == 3 || stress == 4
    if collapse
        if ~exist('Szero','var')
            [Szero] = f_init_stress(Req_zero, Re8, Ca, De, We, Cstar, Pv_star);
        end
    else
        Szero = 0;
    end
elseif stress == 5
    % TODO initial max stress for UCM and Oldroyd-B
    Szero = zeros((Nv - 1)*(spectral == 1) + 2,1);
end

% out parameters

% equation settings
eqns_opts = [radial bubtherm medtherm stress eps3 masstrans];
% solver options
solve_opts = [method spectral divisions Nv Nt Mt Lv Lt];
% dimensionless initial conditions
init_opts = [Rzero Uzero Pb_star P8 T8 Pv_star Req_zero alphax];
% dimensionaless initial stress
init_stress = Szero;
% time span options
tspan_opts = tvector;
% output options
out_opts = [dimensionalout progdisplay];

% physical parameters

% acoustic parameters
acos_opts = [Cstar GAMa kappa nstate hugoniot_s];
% dimensionless waveform parameters
wave_opts = [om ee tw dt mn wave_type wave_poly wave_dpoly];
% dimensionless viscoelastic
sigma_opts = [We Re8 DRe v_a v_nc Ca LAM De JdotA nu_model v_lambda_star zeNO iDRe];
% dimensionless thermal
thermal_opts = [Foh Br alpha_g beta_g alpha_v beta_v chi iota];
% dimensionaless mass transfer
mass_opts = [Fom C0 Rv_star Ra_star L_heat_star mv0 ma0];

end
