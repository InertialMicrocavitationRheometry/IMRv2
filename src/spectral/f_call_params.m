function [vecout]  = f_call_params(varargin)
  % Code to create parameter .mat file for RP_Cav to use 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default options and numerical parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display options
    dimensionalout  = 0;        % output result in dimensional variables
    progdisplay     = 0;        % display progress while code running
    detail          = 2000;    	% number of points in time to store result
    plotresult      = 1;        % generate figure containing results
    radiusonly      = 1;        % only produce R(t) curve
    vitalsreport    = 1;        % display accuracy data
    % output options
    displayonly     = 0;        % do not generate output
    technical       = 0;        % output technical data
    % model for radial bubble dynamics, default is KME in pressure
    % if all values below are zero, Keller-Miksis w/ pressure is used
    rayleighplesset = 1;        % Rayleigh-Plesset equation
    enthalpy        = 0;        % Keller-Miksis enthalpy equations
    gil             = 0;        % Gilmore equation
    % thermal assumptions, default is no thermal assumptions
    polytropic      = 1;        % 0: polytropic assumption, 1: thermal PDE in bubble
    cold            = 1;        % 0: warm fluid, 1: cold fluid assumption
    vapor           = 0;        % 0 : ignore vapor pressure, 1 : vapor pressure
    % mass transfer, default is no mass transfer 
    cgrad           = 0;        % not yet operation leave zero
    % constitutive model
    kelvinVoigt         = 1;        % neo-Hookean
    yangChurch           = 0;        % yangChurch model
    linelas         = 0;        % linear elastic model
    liner           = 0;        % linear Maxwell, Jeffreys, Zener depending on material parameters
    oldb            = 0;        % upper-convected Maxwell, OldRoyd-B depending on material parameters
    ptt             = 0;        % Phan-Thien-Tanner, check material properties
    gies            = 0;        % Giesekus fluid, check material properties
    vmaterial       = 'water';
    % solver options
    method          = 45;       % ode45 setting for the time stepper
    spectral        = 0;        % force spectral collocation solution
    divisions       = 0;        % minimum number of timesteps
    % numerical parameters
    Nt              = 12;       % number of points in bubble, thermal PDE
    Mt              = 12;       % number of points outside of bubble, thermal PDE
    Nv              = 180;      % number of points outside of bubble, viscoelastic stress PDE     
    Lv              = 3;        % characteristic length for grid stretching, leave at 3
    Lt              = 3;        % characteristic length for grid stretching, leave at 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default physical parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R0              = 2E-6;     % initial bubble radius
    U0              = 0;        % initial velocity (m/s)
    Req             = 0.5E-6;     % Equilibrium radius for pre-stress bubble, see Estrada JMPS 2017
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % waveform parameters   %
    %%%%%%%%%%%%%%%%%%%%%%%%% 
    TFin            = 10e-6;     % final time (s)
    pA              = 0*1e5;      % pressure amplitude (Pa)
    omega           = 0*4e6*2*pi; % frequency (rad/s)
    TW              = 0;        % gaussian width (s)
    DT              = 0;        % delay (s)
    mn              = 0;        % power shift for waveform
    wavetype        = 2;        % wave type oscillating bubble, see f_pinfinity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % acoustic parameters, enthalpy values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho8            = 1060;             % far-field density (kg/m^3)   
    GAM             = 3049.13*1e5;      % state equation parameter (Pa)
    nstate          = 7.15;             % state equation parameter
    P8              = 101325;           % far-field pressure (Pa)
    C8              = sqrt(nstate*(P8 + GAM)/rho8); % far-field sound speed (m/s)
    %%%%%%%%%%%%%%%%%%%%%%
    % thermal parameters %
    %%%%%%%%%%%%%%%%%%%%%%
    % gas properties
	kappa           = 1.47;               % Ratio of Specific Heats 
    % medium properties 
    AT              = 5.28e-5;            % (W/m-K^2)Thermal Conductivity coeff
    BT              = 1.165e-2;           % (W/m-K)Thermal Conductivity coeff
	T8              = 300;                % (K) Far field temp. 
    Km              = 0.615;              % (W/m-K)Thermal Conductivity Medium
    Cp              = 4181;               % Specific Heat Medium J/Kg K;
    Dm              = Km / (rho8*Cp) ;    % Thermal Diffusivity m^2/s 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mass diffusion parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D0              = 24.2e-6;            % Diffusion Coeff m^2/s
	L_heat          = 2264.76e3;          % (J/Kg) Latent heat of evaporation
	Ru              = 8.3144598;          % (J/mol-K) Universal Gas Constant % Ru/(18.01528e-3);%Ru/(28.966e-3);     
    Rv              = Ru/(18.01528e-3);   % (J/Kg-K) Gas constant vapor
    Ra              = 438.275;            % (J/Kg-K)Gas constant air    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % constitutive params   %
    %%%%%%%%%%%%%%%%%%%%%%%%%    
    S               = 0.072;              % (N/m) Liquid Surface Tension 
    % non-Newtonian viscosity
    [mu8,Dmu,v_a,v_nc,v_lambda] = f_nonNewtonian_Re(vmaterial); % viscosity
	G               = 1E1;                % (Pa) Medium Shear Modulus 
    lambda1         = 0;  % 0.5e-6;             % relaxation time (s)
    lambda2         = lambda1;            % retardation time (s)
    Pv              = f_pvsat(T8);  
	P0              = Pv*vapor + ...
        (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^(3*kappa); % need to add Pv_sat at room temp  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % overwrite defaults with options and dimensional inputs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check that all inputs are matched
    if mod(nargin,2) == 1
        disp('Error: unmatched inputs'); 
        return; 
    end
    % load inputs
    misscount = 0; p0set = 0; c8set = 0;
    if isempty(varargin) == 1
        disp('Using default settings');
    else
        for n = 1:2:nargin
            if strcmpi(varargin{n},'p0') == 1
                p0set = 1; 
            end
            if strcmpi(varargin{n},'c8') == 1 
                c8set = 1; 
            end
            switch lower(varargin{n})
            % thermal options
            case 'poly',    polytropic = varargin{n+1};% ~= 0;
            case 'cold',    cold = varargin{n+1};% ~= 0;
            case 'vapor',   vapor = varargin{n+1};% ~= 0;
            % equation for radial dynamics
            case 'rpe',     rayleighplesset = varargin{n+1};% ~= 0;
            case 'enth',    enthalpy = varargin{n+1};% ~= 0;
            case 'gil',     gil = varargin{n+1};% ~= 0;
            % constitutive model
            case 'kelvinvoigt', kelvinVoigt = varargin{n+1};% ~= 0;
            case 'yangchurch',   yangChurch = varargin{n+1};% ~= 0;
            case 'linelas', linelas = varargin{n+1};% ~= 0;
            case 'liner',   liner = varargin{n+1};% ~= 0;
            case 'oldb',    oldb = varargin{n+1};% ~= 0;
            case 'ptt',     ptt = varargin{n+1};% ~= 0;
            case 'gies',    gies = varargin{n+1};
            % display options
            case 'dimout',  dimensionalout = varargin{n+1};% ~= 0;
            case 'pdisp',   progdisplay = varargin{n+1};% ~= 0;
            case 'detail',  detail = varargin{n+1};
            case 'plot',    plotresult = varargin{n+1};% ~= 0;
            case 'ronly',   radiusonly = varargin{n+1};% ~= 0;
            case 'vitals',  vitalsreport = varargin{n+1};% ~= 0;
            % output options
            case 'donly',   displayonly = varargin{n+1};% ~= 0;
            case 'tech',    technical = varargin{n+1};% ~= 0;
            % solver options
            case 'method',  method = varargin{n+1};
            case 'spectral',spectral = varargin{n+1};% ~= 0;
            case 'divisions', divisions = varargin{n+1};
            % numerical parameters
            case 'nv',      Nv = varargin{n+1};
            case 'nt',      Nt = varargin{n+1};
            case 'mt',      Mt = varargin{n+1};
            case 'lv',      Lv = varargin{n+1};
            case 'lt',      Lt = varargin{n+1};  
            % waveform parameters
            case 'tfin',    TFin = varargin{n+1};
            case 'pa',      pA = varargin{n+1};
            case 'omega',   omega = varargin{n+1};
            case 'tw',      TW = varargin{n+1};
            case 'dt',      DT = varargin{n+1};
            case 'mn',      mn = varargin{n+1};
            % physical constants
            case 'gam',     GAM = varargin{n+1};
            case 'nstate',  nstate = varargin{n+1};
            case 'p8',      P8 = varargin{n+1};
            case 'rho8',    rho8 = varargin{n+1};
            case 'c8',      C8 = varargin{n+1};
            case 'surf',    S = varargin{n+1};
            case 'kappa',   kappa = varargin{n+1};
            case 't8',      T8 = varargin{n+1};
            case 'pv',      Pv = varargin{n+1};
            % thermal conductivity and diffusivity
            case 'at',      AT = varargin{n+1};
            case 'bt',      BT = varargin{n+1};
            case 'kl',      Km = varargin{n+1};
            case 'dl',      Dm = varargin{n+1};
            % viscoelastic parameters
            case 'mu',      mu8 = varargin{n+1};
            case 'g',       G = varargin{n+1};
            case 'lambda1', lambda1 = varargin{n+1};
            case 'lambda2', lambda2 = varargin{n+1};
            % mass diffusion properties
            case 'dmass',   D0 = varargin{n+1};
            % initial conditions
            case 'r0',      R0 = varargin{n+1};
            %case 'rref',    Rref = varargin{n+1};
            case 'u0',      U0 = varargin{n+1};
            case 'p0',      P0 = varargin{n+1};
            case 'req',     Req = varargin{n+1};
            %case 'a0', a0 = varargin{n+1};
            %case 'b0', b0 = varargin{n+1};
            otherwise, misscount = misscount + 1;
            end
        end
        if p0set == 0
            P0 = P8 + 2*S/Req - Pv*vapor;
        end
        if c8set == 0 
            C8 = sqrt(nstate*(P8 + GAM)/rho8); 
        end
    end
    % check for physical viscoelastic parameters
    if (lambda1 > mu8/G && (yangChurch == 0 && kelvinVoigt == 0 && ...
            linelas == 0)) || abs(gies - 0.25) > 0.25
        disp('Error: nonphysical viscoelastic parameters');
        return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Intermediate calculated variables %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    K8      = AT*T8+BT;                 % far-field thermal conductivity (W/(m K))
    Rnondim = P8/(rho8*T8);             % dimenisonal parameter for gas constants
    Uc      = sqrt(P8/rho8);            % characteristic velocity (m/s)
    theta   = Rv/Ra*(P0-Pv)/Pv;         % mass air / mass vapor 
    C0      = 1/(1+theta);              % initial vapor concentration
    mv0     = Pv*(4/3*pi*R0^3)/Rv/T8;   % mass of vapor
    ma0     = (P0-Pv)*(4/3*pi*R0^3)/Ra/T8;  % mass of non-condensible gas
    Mnondim = rho8*(4/3*pi*R0^3);       %
    
    % Final non-dimensional variables
    t0      = R0/Uc;                    % characteristic time (s) 
	% dimensionless vapor and infinity pressure
    Pv_star = vapor*Pv/P8;
	P0_star = P0/P8;                    % 
    % dimensionless waveform parameters
    tfin    = TFin/t0;                  % simulation time
    om      = omega*t0;                 % non-dimensional frequency
    ee      = pA/P8;                    
    tw      = TW*t0;
    dt      = DT/t0;    
    % acoustic properties
    GAMa    = GAM/P8;                   % bulk liquid stiffness
	Cstar   = C8/Uc;                    % speed of sound
    % thermal properties
    chi     = T8*K8/(P8*R0*Uc);
    iota    = Km/(K8*Lt);
    Foh     = Dm/(Uc*R0); 
    alpha   = AT*T8/K8;              % we do not need beta = BT/K8, we do for diffusion
    beta    = BT/K8;
    Br      = Uc^2/(Cp*T8);      
    % mass diffusion
	Fom     = D0/(Uc*R0);
    mv0     = mv0/Mnondim; 
    ma0     = ma0/Mnondim;     
    Rv_star = Rv/Rnondim;
    Ra_star = Ra/Rnondim; 
    L_heat_star = L_heat/(Uc)^2;
    % viscoelastic properties
    Ca      = P8/G; 
    Re8     = P8*R0/(mu8*Uc);        % this is the Reynolds number
    if Dmu ~= 0
        DRe = P8*R0/(Dmu*Uc);
    else 
        DRe = 0;
    end
    We      = P8*R0/(2*S);
    v_lambda_star = v_lambda/t0;
    LAM     = lambda2/lambda1;
    De      = lambda1*Uc/R0;            % Deborah number
    % dimensionless initial conditions
    Rzero   = 1;
    Req_zero= Req/R0;
    Uzero   = U0/Uc;
    pzero   = P0_star;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overwrite defaults with nondimensional inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin) == 0
    for n = 1:2:nargin
        switch lower(varargin{n})  
            % dimensionless state variables
            case 'cstar',   Cstar = varargin{n+1};
            case 'gama',    GAMa = varargin{n+1};
            % dimensionless waveform parameters
            case 'tfinx',   tfin = varargin{n+1};
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
            case 'alpha',   alpha = varargin{n+1};
            case 'chi',     chi = varargin{n+1};
            case 'iota',    iota = varargin{n+1};
            % dimensionless initial conditions
            case 'rzero',   Rzero = varargin{n+1};
            case 'uzero',   Uzero = varargin{n+1};
            case 'pzero',   pzero = varargin{n+1};
            otherwise
                misscount = misscount + 1;
        end
    end
end
% check that all inputs were accounted for
if misscount ~= nargin/2
    disp(['Error: ' num2str(misscount-nargin/2) ' unrecognized input(s)']);
    return;
end
% check dimensionless viscoelastic parameters
if De > Ca/Re8 && (yangChurch == 0 && kelvinVoigt == 0 && linelas == 0)
    disp('Error: nonphysical viscoelastic parameters');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% final setting adjustments %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rayleighplesset == 1, enthalpy = 0; gil = 0; end
if enthalpy == 1, rayleighplesset = 0; gil = 0; end
if gil == 1, rayleighplesset = 0; enthalpy = 0; end
if polytropic == 1, cold = 0; end
if cold == 1, polytropic = 0; cgrad = 0; end

if kelvinVoigt == 1
    [yangChurch,linelas,liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif yangChurch == 1
    [linelas,liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif linelas == 1
    [liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif liner == 1
    [yangChurch,ptt,gies] = deal(0);
elseif oldb == 1
    [yangChurch,liner,ptt,gies] = deal(0);
    Ca = -1;
elseif ptt == 1
    [yangChurch,liner,gies,LAM] = deal(0);
    Ca = -1; spectral = 1;
elseif gies ~= 0
    [yangChurch,liner,ptt,LAM] = deal(0);
    Ca = -1; spectral = 1;
end

if yangChurch == 1 || kelvinVoigt == 1 || linelas == 1
    JdotA = 4/Re8;
elseif liner == 1 || oldb == 1
    JdotA = 4*LAM/Re8;
else
    JdotA = 0;
end
if spectral == 1
    JdotA = 0; 
end

vecout = {...
      ... % numerical settings 
      polytropic cold cgrad rayleighplesset enthalpy gil ...
      kelvinVoigt yangChurch linelas liner oldb ptt gies ...
      ...% output options
      dimensionalout progdisplay detail plotresult ...
      radiusonly vitalsreport displayonly technical ... 
      ...% solver options
      method spectral divisions... 
      ...% numerical parameters 
      Nv Nt Mt Lv Lt ... 
      ... % physical parameters%
      Cstar GAMa kappa nstate ... % acoustic parameters
      tfin om ee tw dt mn wavetype ... % dimensionless waveform parameters
      We Re8 DRe v_a v_nc Ca LAM De JdotA v_lambda_star ... % dimensionless viscoelastic
      Fom Br alpha beta chi iota ... % dimensionless thermal 
      Foh C0 Rv_star Ra_star L_heat_star mv0 ma0 ... % dimensionaless mass transfer 
      Rzero Uzero pzero P8 T8 Pv_star Req_zero... % dimensionless initial conditions
      };  
end

function [mu8,Dmu,a,nc,lambda] = f_nonNewtonian_Re(vmaterial)
%F_NONNEWTONIAN Outputs the Reynolds number that is dynamically changes
% with the shear rate. Note: Re = P8*R0/(m8*Uc). Units are in Pascal
% seconds.
    a = 0; nc = 0; lambda = 0;
    if strcmp('water',vmaterial)==1
        mu8 = 8.3283e-4;
        muo = 8.3283e-4;
    elseif strcmp('blood_infinity', vmaterial) == 1
        mu8 = 0.00345; 
        muo = mu8; 
    elseif strcmp('blood_zero', vmaterial) == 1
        mu8 = 0.056; 
        muo = mu8; 
    elseif strcmp('blood_combined', vmaterial) == 1
        mu8 = 0.00345; 
        muo = 0.056; 
        nc = 0.384; 
        lambda = 5.61; 
    elseif strcmp('blood_biro', vmaterial) == 1
        mu8 = 0.00345; 
        muo = 0.056; 
        nc = 0.3568; 
        lambda = 2.96;         
    elseif strcmp('blood_merrill', vmaterial) == 1
        mu8 = 0.00345; 
        muo = 0.056; 
        nc = 0.205; 
        lambda = 9.67;         
    elseif strcmp('blood_skalak', vmaterial) == 1
        mu8 = 0.00345; 
        muo = 0.056; 
        nc = 0.218; 
        lambda = 4.58;                 
    elseif strcmp('polystyrene', vmaterial) ==1
        mu8 = 0; 
        muo = 4*10^6;
    elseif strcmp('aluminum soap', vmaterial) == 1
        mu8 = 0.01;
        muo = 89.6; 
    elseif strcmp('p-oxide', vmaterial) == 1
        mu8 = 0; 
        muo = 15.25; 
    elseif strcmp('h-cellulose', vmaterial) == 1
        mu8 = 0; 
        muo = 0.22;
    else
        error('No viscosity model specified in f_call_parameters, exiting');
    end
    Dmu = muo-mu8;
end
