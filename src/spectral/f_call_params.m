function [eqns_opts, solve_opts, init_opts, tspan_opts, out_opts, ...
    acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] ...
    = f_call_params(varargin)
    % Code to create parameter .mat file for RP_Cav to use 

    %*************************************************************************
    % EQUATION OPTIONS
    radial          = 2;        % 1 : RP, 2 : K-M P, 3 : K-M E, 4 : Gil
    bubtherm        = 0;        % 0 : polytropic assumption, 1: thermal PDE in bubble
    medtherm        = 0;        % 0 : cold fluid, 1: warm fluid assumption
    stress          = 1;        % 1 : NHKV, qKV, 2: linear Maxwell, Jeffreys, Zener, 3: UCM or OldB, 4: PTT, 5: Giesekus
    eps3            = 0;        % this value must be (0, 0.5]
    vapor           = 0;        % 0 : ignore vapor pressure, 1 : vapor pressure
    masstrans       = 0;        % mass transfer, default is no mass transfer 
    perturbed       = 0;        % perturbation equation from Jin's paper
    %*************************************************************************
    % SOLVER OPTIONS
    TFin            = 100e-6;     % final time (s)
    TVector         = [0 TFin];
    method          = 23;       % ode45 setting for the time stepper
    spectral        = 0;        % force spectral collocation solution
    divisions       = 0;        % minimum number of timesteps
    Nt              = 12;       % number of points in bubble, thermal PDE
    Mt              = 12;       % number of points outside of bubble, thermal PDE
    Nv              = 150;      % number of points outside of bubble, viscoelastic stress PDE     
    Lv              = 3;        % characteristic length for grid stretching, leave at 3
    Lt              = 3;        % characteristic length for grid stretching, leave at 3
    %*************************************************************************
    % INITIAL CONDITIONS
    R0              = 2.447495043190468e-04;     % initial bubble radius
    U0              = 0;        % initial velocity (m/s)
    Req             = 3.008409399929516e-05;%27E-6;   % Equilibrium radius for pre-stress bubble, see Estrada JMPS 2017    
    a0              = zeros(6,1);
    a0(1)           = 0.00*R0;
    a0(2)           = 0.00*R0;
    a0(3)           = 0.00*R0;
    a0(4)           = 0.00*R0;
    a0(5)           = 0.00*R0;
    a0(6)           = 0.00*R0;
    adot0           = zeros(6,1);
    %*************************************************************************
    % OUPUT OPTIONS*
    dimensionalout  = 0;        % output result in dimensional variables
    progdisplay     = 0;        % display progress while code running
    plotresult      = 0;        % generate figure containing results

    %*************************************************************************
    % ACOUSTIC OPTIONS
    rho8            = 1064;%997;              % far-field density (kg/m^3)   
    GAM             = 3049.13*1e5;      % state equation parameter (Pa)
    nstate          = 7.15;             % state equation parameter
    P8              = 101325;           % far-field pressure (Pa)
    C8              = 1484;%sqrt(nstate*(P8 + GAM)/rho8); % far-field sound speed (m/s)
    %*************************************************************************
    % PRESSURE WAVEFORM OPTIONS
    pA              = 0*1e6;      % pressure amplitude (Pa)
    omega           = 0*4e6*2*pi; % frequency (rad/s)
    TW              = 0;        % gaussian width (s)
    DT              = 0;        % delay (s)
    mn              = 0;        % power shift for waveform
    wavetype        = 2;        % wave type oscillating bubble, see f_pinfinity
    l               = 2:7;     % mode number
    nl              = length(l);
    %*************************************************************************
    % STRESS OPTIONS
    S               = 0.072;%0.056;% 0.072              % (N/m) Liquid Surface Tension 
    vmaterial       = 'water';
    [mu8,Dmu,v_a,v_nc,v_lambda] = f_nonNewtonian_Re(vmaterial); % non-Newtonian viscosity
	G               = 312.5;%1E3;                % (Pa) Medium Shear Modulus 
    lambda1         = 0*0.5e-5;             % relaxation time (s)
    lambda2         = 0;            % retardation time (s)
    mu8             = 0.027606;%0.0246;%1E-3;
    alphax          = 0;        % qKV term
    Ca              = P8/G;           
    S0              = 0;%(3*alphax-1)*(5 - (Req/1)^4 - 4*(Req/1))/(2*Ca) + ...
        % (2*alphax/Ca)*(27/40 + (1/8)*(Req/1)^8 + (1/5)*(Req/1)^5 + (1/2)*(Req/1)^2 - ...
        % 2*1/Req)
    
    %*************************************************************************
    % THERMAL OPTIONS
    % gas properties
	kappa           = 1.4;               % Ratio of Specific Heats 
    % medium properties 
    AT              = 5.28e-5;            % (W/m-K^2)Thermal Conductivity coeff
    BT              = 1.165e-2;           % (W/m-K)Thermal Conductivity coeff
	T8              = 298.15;                % (K) Far field temp. 
    Km              = 0.55;%0.615;              % (W/m-K)Thermal Conductivity Medium
    Cp              = 4181;               % Specific Heat Medium J/Kg K;
    Dm              = Km / (rho8*Cp) ;    % Thermal Diffusivity m^2/s 
	%*************************************************************************
    % MASS TRANSFER OPTIONS
    D0              = 24.2e-6;            % Diffusion Coeff m^2/s
	L_heat          = 2264.76e3;          % (J/Kg) Latent heat of evaporation
	Ru              = 8.3144598;          % (J/mol-K) Universal Gas Constant % Ru/(18.01528e-3);%Ru/(28.966e-3);     
    Rv              = Ru/(18.01528e-3);   % (J/Kg-K) Gas constant vapor
    Ra              = 438.275;            % (J/Kg-K)Gas constant air    
    %*************************************************************************
    % PRESSURE OPTIONS
    Pv              = f_pvsat(T8);  
	% P0              = Pv*vapor + ...% P8 + 2*S/R0;%
    %    (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^(3*kappa); % need to add Pv_sat at room temp  
    P0              = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3))
    %*************************************************************************
    % OVERRIDES defaults with options and dimensional inputs %
    % check that all inputs are matched
    if mod(nargin,2) == 1
        error('INPUT ERROR: unmatched inputs'); 
    end
    % load inputs
    misscount = 0; p0set = 0; c8set = 0; tflag = 0; visflag = 0;

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
            %*****************************************************************
            % EQUATION OPTIONS
            case 'radial',      radial = varargin{n+1};
            case 'bubtherm',    bubtherm = varargin{n+1};
            case 'medtherm',    medtherm = varargin{n+1};
            case 'stress',      stress = varargin{n+1};
            case 'eps3',        eps3 = varargin{n+1};
            case 'vapor',       vapor = varargin{n+1};
            case 'masstrans',   masstrans = varargin{n+1};
            case 'perturbed',   perturbed = varargin{n+1};
            %*****************************************************************
            % SOLVER OPTIONS
            case 'method',  method = varargin{n+1};
            case 'spectral',spectral = varargin{n+1};
            case 'divisions', divisions = varargin{n+1};
            case 'nv',      Nv = varargin{n+1};
            case 'nt',      Nt = varargin{n+1};
            case 'mt',      Mt = varargin{n+1};
            case 'lv',      Lv = varargin{n+1};
            case 'lt',      Lt = varargin{n+1};  
            case 'tfin',    TFin = varargin{n+1};  
                tflag = tflag + 1; %TVector = 0;
            %*****************************************************************
            % INITIAL CONDITIONS
            case 'r0',      R0 = varargin{n+1};
                P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3))
            case 'u0',      U0 = varargin{n+1};
            case 'req',     Req = varargin{n+1};
                P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3))
            case 'a0',      a0 = varargin{n+1};
            case 'adot0',   adot0 = varargin{n+1}; 
            case 's0',      S0 = varargin{n+1};
            %*****************************************************************
            % OUPUT OPTIONS
            case 'dimout',  dimensionalout = varargin{n+1};
            case 'pdisp',   progdisplay = varargin{n+1};
            case 'plot',    plotresult = varargin{n+1};

            %*****************************************************************
            % ACOUSTIC OPTIONS
            case 'rho8',    rho8 = varargin{n+1};
            case 'gam',     GAM = varargin{n+1};
            case 'nstate',  nstate = varargin{n+1};
            case 'p8',      P8 = varargin{n+1};
            case 'c8',      C8 = varargin{n+1};
            %*****************************************************************
            % PRESSURE WAVEFORM OPTIONS
            case 'pa',      pA = varargin{n+1};
            case 'omega',   omega = varargin{n+1};
            case 'tw',      TW = varargin{n+1};
            case 'dt',      DT = varargin{n+1};
            case 'mn',      mn = varargin{n+1};
            %****************************************************************
            % STRESS OPTIONS
            case 'mu',      mu8 = varargin{n+1};  
                visflag = visflag + 1;
            case 'g',       G = varargin{n+1};
            case 'lambda1', lambda1 = varargin{n+1};
            case 'lambda2', lambda2 = varargin{n+1};
            case 'alphax', alphax = varargin{n+1};
            case 'surft',    S = varargin{n+1};
                P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
            case 'vmaterial', vmaterial = varargin{n+1}; 
                visflag = visflag + 1;
               [mu8,Dmu,v_a,v_nc,v_lambda] = f_nonNewtonian_Re(vmaterial); % non-Newtonian viscosity    
            %*****************************************************************
            % THERMAL OPTIONS            
            case 't8',      T8 = varargin{n+1};
            case 'kappa',   kappa = varargin{n+1};                
            case 'at',      AT = varargin{n+1};
            case 'bt',      BT = varargin{n+1};
            case 'kl',      Km = varargin{n+1};
            case 'dl',      Dm = varargin{n+1};
	        %*****************************************************************
            % MASS TRANSFER OPTIONS
            case 'dmass',   D0 = varargin{n+1};
            case 'lheat',   L_heat = varargin{n+1};
            case 'rv',      Rv = varargin{n+1};
            case 'ra',      Ra = varargin{n+1};
            %*****************************************************************
            % PRESSURE OPTIONS         
            case 'pv',      Pv = varargin{n+1}; 
                P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
            case 'p0',      P0 = varargin{n+1};
            case 'tvector', TVector = varargin{n+1}; 
                tflag = tflag + 1; TFin = 0;
            otherwise, misscount = misscount + 1;
            end
        end
        if p0set == 0
          % need to add Pv_sat at room temp  
          P0 = (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^(3*kappa);
        end
        if c8set == 0 
            C8 = sqrt(nstate*(P8 + GAM)/rho8); 
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
        if check || stress > 5 || stress <= 0
            error('INPUT ERROR: stress must be 1, 2, 3, 4, or 5');
        end        
        check = 1-isnumeric(vapor);
        if check || vapor ~= 0 && vapor ~= 1
            error('INPUT ERROR: vapor must be 0 or 1');
        end
        check = 1-isnumeric(radial);
        if check || masstrans ~= 0 && masstrans ~= 1
            error('INPUT ERROR: masstrans must be 0 or 1');
        end        
        if (tflag > 1)
            error('INPUT ERROR: Only tvector or tfin can be specified, not both');
        end
        if (visflag > 1)
            error('INPUT ERROR: Only vmaterial or mu8 can be specified, not both');
        end        
        if TVector == 0
            TVector = [0 TFin];
        end
    end
    % check for physical viscoelastic parameters
    if (lambda1 > mu8/G && (stress == 1)) || abs(eps3 - 0.25) > 0.25
        error('INPUT ERROR: Unphysical viscoelastic parameters');
    end

    %*************************************************************************
    % Intermediate calculated variables 
    K8      = AT*T8+BT;                 % far-field thermal conductivity (W/(m K))
    Rnondim = P8/(rho8*T8);             % dimenisonal parameter for gas constants
    Uc      = sqrt(P8/rho8);            % characteristic velocity (m/s)
    theta   = Rv/Ra*(P0-Pv)/Pv;         % mass air / mass vapor 
    C0      = 1/(1+theta);              % initial vapor concentration
    mv0     = Pv*vapor*(4/3*pi*R0^3)/Rv/T8;   % mass of vapor
    ma0     = (P0-Pv*vapor)*(4/3*pi*R0^3)/Ra/T8;  % mass of non-condensible gas
    Mnondim = rho8*(4/3*pi*R0^3);       %
    
    % Final non-dimensional variables
    Pref    = P8;
    t0      = R0/Uc;                    % characteristic time (s) 
	% dimensionless vapor and infinity pressure
    Pv_star = vapor*Pv/Pref;
	P0_star = P0/Pref;                    % 
    % dimensionless waveform parameters
    tvector = TVector./t0;
    om      = omega*t0;                 % non-dimensional frequency
    ee      = pA/Pref;                    
    tw      = TW*t0;
    dt      = DT/t0;    
    % acoustic properties
    GAMa    = GAM/P8;                   % bulk liquid stiffness
	Cstar   = C8/Uc;                    % speed of sound
    % thermal properties
    chi     = T8*K8/(P8*R0*Uc);% T8*K8/(P0*R0*Uc); 
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
    Ca      = Pref/G; 
    Re8     = Pref*R0/(mu8*Uc);        % this is the Reynolds number
    if Dmu ~= 0
        DRe = Pref*R0/(Dmu*Uc);
    else 
        DRe = 0;
    end
    We      = Pref*R0/(2*S);
    v_lambda_star = v_lambda/t0;
    LAM     = lambda2/lambda1;
    De      = lambda1*Uc/R0;            % Deborah number
    % dimensionless initial conditions
    Rzero   = 1;
    Req_zero= Req/R0
    Uzero   = U0/Uc;
    pzero   = P0_star;
    azero   = a0./R0;
    adot_zero=adot0./Uc;

%*************************************************************************
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
    error(['INPUT ERROR: ' num2str(misscount-nargin/2) ' unrecognized input(s)']);
end
% check dimensionless viscoelastic parameters
% if De > Ca/Re8 && (yangChurch == 0 && kelvinVoigt == 0 && linelas == 0)
%     error('INPUT ERROR: nonphysical viscoelastic parameters');
% end

%*************************************************************************
% final setting adjustments 
if bubtherm == 0, medtherm = 0; end
if medtherm == 1, bubtherm = 1; masstrans = 0; end

% 1 : N-H, 2: linear Maxwell, Jeffreys, Zener, 3: UCM or OldB, 4: PTT, 5: Giesekus
if stress == 1
    spectral = 0;
elseif stress == 2
    spectral = 0;
elseif stress == 3
    Ca = -1; 
elseif stress == 4
    Ca = -1; spectral = 1;
elseif stress == 5
    Ca = -1; spectral = 1;
end

if stress == 1 
    JdotA = 4/Re8;
elseif stress == 2 || stress == 3 
    JdotA = 4*LAM/Re8;
else
    JdotA = 0;
end
if spectral == 1
    JdotA = 0; 
end

%*************************************************************************
% equation settings 
eqns_opts = [radial bubtherm medtherm stress eps3 vapor masstrans perturbed nl];
% solver options
solve_opts = [method spectral divisions Nv Nt Mt Lv Lt];
% dimensionless initial conditions
init_opts = [Rzero Uzero pzero P8 T8 Pv_star Req_zero S0 alphax azero' adot_zero'];
% time span options
tspan_opts = tvector;
% output options
out_opts = [dimensionalout progdisplay plotresult];

% physical parameters%%%%
% acoustic parameters
acos_opts = [Cstar GAMa kappa nstate];
% dimensionless waveform parameters
wave_opts = [om ee tw dt mn wavetype l];
% dimensionless viscoelastic
sigma_opts = [We Re8 DRe v_a v_nc Ca LAM De JdotA v_lambda_star];
% dimensionless thermal 
thermal_opts = [Fom Br alpha beta chi iota];
% dimensionaless mass transfer 
mass_opts = [Foh C0 Rv_star Ra_star L_heat_star mv0 ma0];

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
        error('INPUT ERROR: No viscosity model specified in f_call_parameters, exiting');
    end
    Dmu = muo-mu8;
end