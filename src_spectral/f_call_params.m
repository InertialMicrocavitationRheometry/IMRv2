function [varargout]  = f_call_params(varargin)
  % Code to create parameter .mat file for RP_Cav to use 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default options and numerical parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display options
    dimensionalout = 0; % output result in dimensional variables
    progdisplay = 1; % display progress while code running
    detail = 2000; % number of points in time to store result
    plotresult = 1; % generate figure containing results
    radiusonly = 0; % only produce R(t) curve
    vitalsreport = 1; % display accuracy data
    % output options
    displayonly = 0; % do not generate output
    technical = 0; % output technical data
    % model for radial bubble dynamics, default is KME in pressure
    rayleighplesset = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    enthalpy = 0;
    % thermal assumptions, default is no thermal assumptions
    polytropic = 1;
    cold = 0;
    vapor = 0; % ignore vapor pressure
    % constitutive model, default is UCM with linear elasticity
    neoHook = 0;
    voigt = 1;
    linelas = 0;
    liner = 0;
    oldb = 0;
    ptt = 0;
    gies = 0;
    % solver options
    method = '45'; % ode45
    spectral = 0; % force spectral collocation solution
    divisions = 0; % minimum number of timesteps
    % numerical parameters
    Nv = 20; 
    Nt = 10; 
    Mt = 10; 
    Lv = 3; 
    Lt = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default physical parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Parameters: 
    A     = 5.28e-5;            % (W/m-K^2)Thermal Conductivity coeff
    B     = 1.165e-2;           % (W/m-K)Thermal Conductivity coeff
    D0    = 24.2e-6;            % Diffusion Coeff m^2/s
    k     = 1.47;               % Ratio of Specific Heats 
    S     = 0.056;              % (N/m) Liquid Surface Tension 
    G     = 0*(1000)*1E3;       % (Pa) Medium Shear Modulus 
    T_inf = 300;                % (K) Far field temp. 
    P_inf = 101325;             % (Pa) Atmospheric Pressure 
    rho   = 999;                %1050.761; % (Kg/m^3) Liquid Density
    Km    = 0.615;              % (W/m-K)Thermal Conductivity Medium
    Cp    = 3.61e3;             % Specific Heat Medium J/Kg K;
    Dm    = Km /(rho*Cp) ;      % Thermal Diffusivity m^2/s 
    Ru    = 8.3144598;          % (J/mol-K) Universal Gas Constant
    Rv    = Ru/(18.01528e-3);   % (J/Kg-K) Gas constant vapor
    Ra    = 438.275;            % (J/Kg-K)Gas constant air
    % Ru/(18.01528e-3);%Ru/(28.966e-3); 
    L     = 1;                  % Strech variable to map domain outside the bubble
    L_heat= 2264.76e3;          % (J/Kg) Latent heat of evaporation 
    C     = 1510;               % sound speed (m/s)
    [mu8,Dmu,v_a,v_nc,v_lambda] = f_nonNewtonian_Re(vmaterial); % viscosity
    
    % Intermediate calculated variables
    K_infy  = A*T_inf+B; 
    Rnondim = P_inf/(rho*T_inf);
    Uc      = sqrt(P_inf/rho);
    Pv      = f_pvsat(T_inf );
    P0      = P_inf + 2*S/R0 ;      % need to add Pv_sat at room temp 
    thetha  = Rv/Ra*(P0-Pv)/Pv;     % mass air / mass vapor 
    C0      = 1/(1+thetha); 
    mv0     = Pv*(4/3*pi*R0^3)/Rv/T_inf;
    ma0     = (P0-Pv)*(4/3*pi*R0^3)/Ra/T_inf;
    Mnondim = rho*(4/3*pi*R0^3);
    
     
    
    % waveform parameters
    TFin = 0.6e-6; % final time (s)
    pA = 2e6; % pressure amplitude (Pa)
    omega = 4e6*2*pi; % frequency (rad/s)
    TW = 0; % gaussian width (s)
    DT = 0; % delay (s)
    
    % physical constants
    GAM = 3049.13*1e5; % state equation parameter (Pa)
    nstate = 7.15; % state equation parameter
    p8 = 101325; % far-field pressure (Pa)
    rho8 = 1060; % far-field density (kg/m^3)
    c8 = sqrt(nstate*(p8 + GAM)/rho8); % far-field sound speed (m/s)
    kappa = 1.4; % ratio of specific heats
    S = 0.072; % surface tension (N/m)
    T8 = 293.15; % far-field temperature (K)
    cpval = 4181; % specific heat of medium (J/(kg K))
    pV = 2300;%2339.262; % vapor pressure of water (Pa)

    % thermal conductivity and diffusivity
    AT = 5.28e-5; % gas thermal conductivity slope parameter (W/(m K^2))
    BT = 1.165e-2; % gas thermal conductivity shift parameter (W/(m K))
    KL = 0.55; % medium thermal conductivity (W/(m K))
    DL = 1.41e-7; % medium thermal diffusivity
    
    % viscoelastic parameters
    mu = 0.03; % viscosity (Pa s)
    G = 20e3; % shear modulus (Pa) 
    lambda1 = 0.5e-6; % relaxation time (s)
    lambda2 = lambda1; % retardation time (s)

    % initial conditions
    R0 = 1e-6; % initial radius (m)
    U0 = 0; % initial velocity (m/s)
    p0 = p8 + 2*S/R0 - pV*vapor; % initial pressure (Pa)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % overwrite defaults with options and dimensional inputs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check that all inputs are matched
    if mod(nargin,2) == 1, disp('Error: unmatched inputs'); return; end
    % load inputs
    misscount = 0; p0set = 0; c8set = 0;
    if isempty(varargin) == 1
        disp('Using default settings');
    else
        for n = 1:2:nargin
            if strcmpi(varargin{n},'p0') == 1, p0set = 1; end
            if strcmpi(varargin{n},'c8') == 1, c8set = 1; end
            switch lower(varargin{n})
            % thermal options
            case 'poly', polytropic = varargin{n+1} ~= 0;
            case 'cold', cold = varargin{n+1} ~= 0;
            case 'vapor', vapor = varargin{n+1} ~= 0;
            % equation for radial dynamics
            case 'rpe', rayleighplesset = varargin{n+1} ~= 0;
            case 'enth', enthalpy = varargin{n+1} ~= 0;
            case 'gil', gil = varargin{n+1} ~= 0;
            % constitutive model
            case 'neohook', neoHook = varargin{n+1} ~= 0;
            case 'voigt', voigt = varargin{n+1} ~= 0;
            case 'linelas', linelas = varargin{n+1} ~= 0;
            case 'liner', liner = varargin{n+1} ~= 0;
            case 'oldb', oldb = varargin{n+1} ~= 0;
            case 'ptt', ptt = varargin{n+1} ~= 0;
            case 'gies', gies = varargin{n+1};
            % display options
            case 'dimout', dimensionalout = varargin{n+1} ~= 0;
            case 'pdisp', progdisplay = varargin{n+1} ~= 0;
            case 'detail', detail = varargin{n+1};
            case 'plot', plotresult = varargin{n+1} ~= 0;
            case 'ronly', radiusonly = varargin{n+1} ~= 0;
            case 'vitals', vitalsreport = varargin{n+1} ~= 0;
            % output options
            case 'donly', displayonly = varargin{n+1} ~= 0;
            case 'tech', technical = varargin{n+1} ~= 0;
            % solver options
            case 'method', method = varargin{n+1};
            case 'spectral', spectral = varargin{n+1} ~= 0;
            case 'divisions', divisions = varargin{n+1};
            % numerical parameters
            case 'nv', Nv = varargin{n+1};
            case 'nt', Nt = varargin{n+1};
            case 'mt', Mt = varargin{n+1};
            case 'lv', Lv = varargin{n+1};
            case 'lt', Lt = varargin{n+1};  
            % waveform parameters
            case 'tfin', TFin = varargin{n+1};
            case 'pa', pA = varargin{n+1};
            case 'omega', omega = varargin{n+1};
            case 'tw', TW = varargin{n+1};
            case 'dt', DT = varargin{n+1};
            case 'mn', mn = varargin{n+1};
            % physical constants
            case 'gam', GAM = varargin{n+1};
            case 'nstate', nstate = varargin{n+1};
            case 'p8', p8 = varargin{n+1};
            case 'rho8', rho8 = varargin{n+1};
            case 'c8', c8 = varargin{n+1};
            case 'kappa', kappa = varargin{n+1};
            case 'surf', S = varargin{n+1};
            case 't8', T8 = varargin{n+1};
            case 'pv', pV = varargin{n+1};
            % thermal conductivity and diffusivity
            case 'at', AT = varargin{n+1};
            case 'bt', BT = varargin{n+1};
            case 'kl', KL = varargin{n+1};
            case 'dl', DL = varargin{n+1};
            % viscoelastic parameters
            case 'mu', mu = varargin{n+1};
            case 'g', G = varargin{n+1};
            case 'lambda1', lambda1 = varargin{n+1};
            case 'lambda2', lambda2 = varargin{n+1};
            % initial conditions
            case 'r0', R0 = varargin{n+1};
            case 'rref', Rref = varargin{n+1};
            case 'u0', U0 = varargin{n+1};
            case 'p0', p0 = varargin{n+1};
            %case 'a0', a0 = varargin{n+1};
            %case 'b0', b0 = varargin{n+1};
            otherwise, misscount = misscount + 1;
            end
        end
        if p0set == 0, p0 = p8 + 2*S/R0 - pV*vapor; end
        if c8set == 0, c8 = sqrt(nstate*(p8 + GAM)/rho8); end
    end

    % check for physical viscoelastic parameters
    if (lambda1 > mu/G && (voigt == 0 && neoHook == 0 && ...
            linelas == 0)) || abs(gies - 0.25) > 0.25
        disp('Error: nonphysical viscoelastic parameters');
        return;
    end
    
    % Final non-dimensional variables
    chi = T_inf*K_infy/(P_inf*R0*Uc);
    fom = D0/(Uc*R0);
    foh = Dm/(Uc*R0); 
    Ca = P_inf/G; 
    Re8 = P_inf*R0/(mu8*Uc);
    if Dmu ~= 0
        DRe = P_inf*R0/(Dmu*Uc);
    else 
        DRe = 0;
    end
    We = P_inf*R0/(2*S);
    Br = Uc^2/(Cp*T_inf);   
    A_star = A*T_inf /  K_infy;
    B_star = B / K_infy; 
    Rv_star = Rv/Rnondim;
    Ra_star = Ra/Rnondim;
    P0_star = P0/P_inf;
    t0 = R0/Uc;
    L_heat_star = L_heat/(Uc)^2;
    Km_star = Km/K_infy; 
    C_star = C/Uc; 
    mv0 = mv0/Mnondim; ma0 = ma0/Mnondim;     
    v_lambda_star = v_lambda/t0;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nondimensionalize problem %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derived constants
uc = sqrt(p8/rho8); % characteristic velocity (m/s)
tc = Rref/uc; % characteristic time (s)
K8 = AT*T8 + BT; % far-field thermal conductivity (W/(m K))
% dimensionless state variables
C = c8/uc;
GAMa = GAM/p8;
% dimensionless waveform parameters
tfin = TFin/tc;
om = omega*tc;
ee = pA/p8;
tw = TW*tc;
dt = DT/tc;
% dimensionless vapor and infinity pressure
pVap = pV/p8*vapor;
p8bar = p8/p8;
p0star = p0/p8;
% dimensionless numbers
We = (2*S)/(p8*Rref); % Weber number
Re = p8*Rref/(mu*uc); % Reynolds number
Ca = p8/G; % Cauchy number
LAM = lambda2/lambda1;
De = lambda1*uc/Rref; % Deborah number
Fo = DL/(uc*Rref); % Fourier number
Br = uc^2/(T8*cpval); % Brinkman number
% dimensionless thermal quantities
alpha = AT*T8/K8;
chi = T8*K8/(p8*Rref*uc);
iota = KL/(K8*Lt);
% dimensionless initial conditions
Rzero = R0/Rref;
Uzero = U0/uc;
pzero = (p0star)*(R0/Rref)*(3*kappa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overwrite defaults with nondimensional inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin) == 0
    for n = 1:2:nargin
        switch lower(varargin{n})  
            % dimensionless state variables
            case 'cx', C = varargin{n+1};
            case 'gama', GAMa = varargin{n+1};
            % dimensionless waveform parameters
            case 'tfinx', tfin = varargin{n+1};
            case 'om', om = varargin{n+1};
            case 'ee', ee = varargin{n+1};
            case 'twx', tw = varargin{n+1};
            case 'dtx', dt = varargin{n+1}; 
            % dimensionless numbers
            case 'we', We = varargin{n+1};
            case 're', Re = varargin{n+1};
            case 'ca', Ca = varargin{n+1};
            case 'lam', LAM = varargin{n+1};
            case 'de', De = varargin{n+1};
            case 'fo', Fo = varargin{n+1};
            % dimensionless thermal quantities
            case 'alpha', alpha = varargin{n+1};
            case 'chi', chi = varargin{n+1};
            case 'iota', iota = varargin{n+1};
            % dimensionless initial conditions
            case 'rzero', Rzero = varargin{n+1};
            case 'uzero', Uzero = varargin{n+1};
            case 'pzero', pzero = varargin{n+1};
            otherwise, misscount = misscount + 1;
        end
    end
end
% check that all inputs were accounted for
if misscount ~= nargin/2
    disp(['Error: ' num2str(misscount-nargin/2) ' unrecognized input(s)']);
    return;
end
% check dimensionless viscoelastic parameters
if De > Ca/Re && (voigt == 0 && neoHook == 0 && linelas == 0)
    disp('Error: nonphysical viscoelastic parameters');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% final setting adjustments %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rayleighplesset == 1, enthalpy = 0; end
if enthalpy == 1, rayleighplesset = 0; end
if polytropic == 1, cold = 0; end
if cold == 1, polytropic = 0; end

if neoHook == 1
    [voigt,linelas,liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif voigt == 1
    [linelas,liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif linelas == 1
    [liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif liner == 1
    [voigt,ptt,gies] = deal(0);
elseif oldb == 1
    [voigt,liner,ptt,gies] = deal(0);
    Ca = Inf;
elseif ptt == 1
    [voigt,liner,gies,LAM] = deal(0);
    Ca = Inf; spectral = 1;
elseif gies ~= 0
    [voigt,liner,ptt,LAM] = deal(0);
    Ca = Inf; spectral = 1;
end

if voigt == 1 || neoHook == 1 || linelas == 1
    JdotA = 4/Re;
elseif liner == 1 || oldb == 1
    JdotA = 4*LAM/Re;
else
    JdotA = 0;
end
if spectral == 1, JdotA = 0; end
    

  P = [k chi fom foh Ca Re8 We Br A_star...
         B_star Rv_star Ra_star P0_star t0 C0 L L_heat_star Km_star ...
         P_inf  T_inf C_star mv0  ma0 DRe v_a v_nc v_lambda_star];
     
  varargout = {polytropic cold rayleighplesset enthalpy gil ...
      neoHook voigt linelas liner oldb ptt gies ...
      dimensionalout progdisplay detail plotresult radiusonly vitalsreport...
      displayonly technical... % output options
      method spectral divisions... % solver options
      Nv Nt Mt Lv Lt... % numerical parameters 
      C GAMa... % state variables
      tfin om ee tw dt ... % dimensionless waveform parameters
      We Re Ca LAM De Fo ... % dimensionless numbers
      alpha chi iota... % dimensionless thermal quantities
      Rzero Uzero pzero ... % dimensionless initial conditions
      };  
end

function [mu8,Dmu,a,nc,lambda] = f_nonNewtonian_Re(vmaterial)
%F_NONNEWTONIAN Outputs the Reynolds number that is dynamically changes
% with the shear rate. Note: Re = P_inf*R0/(m8*Uc). Units are in Pascal
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