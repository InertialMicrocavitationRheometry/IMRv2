% file f_call_params.m
% brief contains function f_call_params

% brief This function sets up the simulation environment. The function can
% receive an input file. The input file must be in the formatted similar to
% the default_case.m file. Otherwise, if no input file is provided, the
% code will read in the default case. The default case can be altered based
% on the input arguments when invoking either the finite difference or
% spectral module.
function [eqns_opts, solve_opts, init_opts, tspan_opts, out_opts, ...
        acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] ...
    = f_call_params(varargin)
    
    disp('--- Inertial Microcavitation Rheometry Forward Solver ---');
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
    c8set = 0;
    tflag = 0;
    visflag = 0;
    for n = 1:2:nargin
        if strcmpi(varargin{n},'p0') == 1
            p0set = 1;
        end
        if strcmpi(varargin{n},'c8') == 1
            c8set = 1;
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
            case 'method',  method = varargin{n+1};
            case 'spectral',spectral = varargin{n+1};
            case 'divisions', divisions = varargin{n+1};
            case 'nv',      Nv = varargin{n+1};
            case 'nt',      Nt = varargin{n+1};
            case 'mt',      Mt = varargin{n+1};
            case 'lv',      Lv = varargin{n+1};
            case 'lt',      Lt = varargin{n+1};
            case 'tfin',    TFin = varargin{n+1};
            tflag = tflag + 1;
            %TVector = 0;
            
            % initial options
            case 'r0',      R0 = varargin{n+1};
            P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
            case 'u0',      U0 = varargin{n+1};
            case 'req',     Req = varargin{n+1};
            P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
            
            % output options
            case 'dimout',     dimensionalout = varargin{n+1};
            case 'progdisplay',   progdisplay = varargin{n+1};
            
            % acoustic options
            case 'rho8',    rho8 = varargin{n+1};
            case 'gam',     GAM = varargin{n+1};
            case 'nstate',  nstate = varargin{n+1};
            case 'p8',      P8 = varargin{n+1};
            case 'c8',      C8 = varargin{n+1};
            
            % pressure waveform options
            case 'pa',      pA = varargin{n+1};
            case 'omega',   omega = varargin{n+1};
            case 'tw',      TW = varargin{n+1};
            case 'dt',      DT = varargin{n+1};
            case 'mn',      mn = varargin{n+1};
            
            % stress options
            case 'mu',      mu8 = varargin{n+1};
            visflag = visflag + 1;
            
            case 'g',       G = varargin{n+1};
            case 'lambda1', lambda1 = varargin{n+1};
            case 'lambda2', lambda2 = varargin{n+1};
            case 'alphax',  alphax = varargin{n+1};
            case 'surft',    S = varargin{n+1};
            % P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
            case 'vmaterial', vmaterial = varargin{n+1};
            visflag = visflag + 1;
            [mu8,Dmu,v_a,v_nc,v_lambda,vmat] = f_nonNewtonian_Re(vmaterial); % non-Newtonian viscosity
            
            % thermal options
            case 't8',      T8 = varargin{n+1};
            case 'kappa',   kappa = varargin{n+1};
            case 'at',      AT = varargin{n+1};
            case 'bt',      BT = varargin{n+1};
            case 'kl',      Km = varargin{n+1};
            case 'dl',      Dm = varargin{n+1};
            
            % mass transfer options
            case 'dmass',   D0 = varargin{n+1};
            case 'lheat',   L_heat = varargin{n+1};
            case 'rv',      Rv = varargin{n+1};
            case 'ra',      Ra = varargin{n+1};
            
            % pressure options
            case 'pv',      Pv = varargin{n+1};
            P0 = (P8 + 2*S/Req - Pv*vapor)*((Req/R0)^(3));
            case 'p0',      P0 = varargin{n+1};
            case 'tvector', TVector = varargin{n+1};
            tflag = tflag + 1;
            TFin = 0;
            otherwise, misscount = misscount + 1;
            
        end
    end

    if p0set == 0
        % need to add Pv_sat at room temp
        P0 = (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^(3);
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
    
    % check for physical viscoelastic parameters
    if (lambda1 > mu8/G && (stress == 3 || stress == 4)) || abs(eps3 - 0.25) > 0.25
        error('INPUT ERROR: Unphysical viscoelastic parameters');
    end
    
    % intermediate calculated variables
    K8      = AT*T8+BT;                 % far-field thermal conductivity (W/(m K))
    Rnondim = P8/(rho8*T8);             % dimensional parameter for gas constants
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
    P0_star = P0/Pref;
    
    if Pv_star > P0_star
        error('vapor pressure can not be higher than the total pressure');
    end

    % TODO
    % When we assume water vapor undergoes infinitely fast mass diffusion
    % the vapor pressure is constant and P is the pressure of
    % non-condesible gas
    %P0_star = P0_star - (1-masstrans)*f_pvsat(T_inf)/P_inf;
    
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
    chi     = T8*K8/(P8*R0*Uc);
    % T8*K8/(P0*R0*Uc);
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
    log_We = log(0.5) + log(Pref) + log(R0) - log(S);
    We = exp(log_We);
    v_lambda_star = v_lambda/t0;
    LAM     = lambda2/lambda1;
    De      = lambda1*Uc/R0;            % Deborah number
    % dimensionless initial conditions
    Rzero   = 1;
    Req_zero= Req/R0;
    Uzero   = U0/Uc;
    pzero   = P0_star;
    
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
    
    % final setting adjustments
    if bubtherm == 0
        medtherm = 0;
    end
    if medtherm == 1
        bubtherm = 1;
    end
    
    % 1 : N-H, 2: qN-H, 3: linear Maxwell, Jeffreys, Zener, 5: UCM or OldB, 6: PTT, 7: Giesekus
    if stress == 1 || stress == 2 || stress == 3 || stress == 4
        spectral = 0;
    elseif stress == 5
        Ca = -1;
    elseif stress == 6
        Ca = -1;
        spectral = 1;
    elseif stress == 7
        Ca = -1;
        spectral = 1;
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
    if stress == 0 || stress == 1 || stress == 2 || stress == 3 || stress == 4
        zeNO = 0;
    else
        zeNO = 1;
    end
    
    % TODO
    % need to modify initial conditions for the Out-of-Equilibrium Rayleigh
    % collapse:
    % if  (Pext_type == 'IC')
    %     Pv = f_pvsat(1*T_inf)/P_inf;
    %     P0_star = Pext_Amp_Freq(1)/P_inf + Cgrad*f_pvsat(1*T_inf)/P_inf;
    %     % Need to recalculate initial concentration
    %     theta = Rv_star/Ra_star*(P0_star-Pv)/Pv; % masp air / mass vapor
    %     C0 = 1/(1+theta);
    %     Ma0 = (P0_star-Pv)/Ra_star;
    %     % Need to calculate the equilibrium radii for initial
    %     % stress state:
    %     [REq] = f_calc_Req(R0, Tgrad ,Cgrad,Pext_Amp_Freq(1),vmaterial);
    %     %REq = 1;
    %     C0 = C0*ones(1,NT);
    %     U0_star = -(1-P0_star)/(C_star); % Initial velocity
    %     %Plesset & Prosperetti, ARFM 1977, p166
    %     else
    %     U0_star = 0;
    % end
    %
    % if  (Pext_type == 'RC')
    %     U0_star = -1*(Pext_Amp_Freq(1)/P_inf)/(C_star); % Initial velocity
    %      %Plesset & Prosperetti, ARFM 1977, p166
    % end
    
    % equation settings
    eqns_opts = [radial bubtherm medtherm stress eps3 vapor masstrans];
    % solver options
    solve_opts = [method spectral divisions Nv Nt Mt Lv Lt];
    % dimensionless initial conditions
    init_opts = [Rzero Uzero pzero P8 T8 Pv_star Req_zero alphax];
    % time span options
    tspan_opts = tvector;
    % output options
    out_opts = [dimensionalout progdisplay];
    
    % physical parameters
    
    % acoustic parameters
    acos_opts = [Cstar GAMa kappa nstate];
    % dimensionless waveform parameters
    wave_opts = [om ee tw dt mn wave_type];
    % dimensionless viscoelastic
    sigma_opts = [We Re8 DRe v_a v_nc Ca LAM De JdotA vmat v_lambda_star zeNO];
    % dimensionless thermal
    thermal_opts = [Foh Br alpha beta chi iota];
    % dimensionaless mass transfer
    mass_opts = [Fom C0 Rv_star Ra_star L_heat_star mv0 ma0]; 
end