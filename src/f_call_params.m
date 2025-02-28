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
            case 'wave_type', wave_type = varargin{n+1};
            
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
            tempset = 1;
            case 'kappa',   kappa = varargin{n+1};
            case 'at',      AT = varargin{n+1};
            case 'bt',      BT = varargin{n+1};
            case 'km',      Km = varargin{n+1};
            case 'dm',      Dm = varargin{n+1}; 
                            dmset = 1;
            
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
    
    if tempset == 1
        % recalculating the vapor pressure
        Pv = vapor*f_pvsat(T8);
        P0 = (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^3;
    end
    if p0set == 0
        % need to add Pv_sat at room temp
        P0 = (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^3;
    end
    if dmset == 0
        % (m^2/s) thermal diffusivity
        Dm = Km / (rho8*Cp);
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
    check = 1-isnumeric(wave_type);
    if check || wave_type > 5 || wave_type <= -3
        error('INPUT ERROR: wavetype must be between -2 to 5');
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
    
    % far-field thermal conductivity (W/(m K))
    K8      = AT*T8+BT;                 
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
    alpha   = AT*T8/K8;              
    beta    = BT/K8;
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
    Req_zero= Req/R0;
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
                case 'alpha',   alpha = varargin{n+1};
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
            % Keller-Miksis equation
            % if linkv==1 || neoHook==1 || Yeoh==1
                % JdotA = 4/Re8;
            % elseif sls==1 || nhzen==1 || fdkv==1 || zzzen==1 || fdmax==1  % ZZ - Not exactly true for FDKV, but I also don't think this matters ...
                % JdotA = 0;
            % elseif nhkv_pld==1
            %     %JdotA = 4/Re8*(2^alpha+1)/3*(abs(Rdot)/R)^(alpha-1);
            %     JdotA = 4/Re8/3*(2^alpha+1)*sign(Rdot)*(abs(Rdot)/R)^(alpha)*R^2/Rdot^2;
            %     if isnan(SdotA)
            %         SdotA=4/Re8;
            %     end
            % end

    if stress == 0 || stress == 1 || stress == 2 || stress == 3 || stress == 4
        zeNO = 0;
    else
        zeNO = 1;
    end
    
    % modify initial conditions for the Out-of-Equilibrium Rayleigh collapse:
    % if  (wave_type == -1)
    %     % mass air / mass vapor
    %     theta = Rv_star/Ra_star*P0/Pv; 
    %     C0 = 1/(1+theta);
    %     ma0 = P0/Ra_star;
    %     % calculating the equilibrium radii for initial
    %     [REq] = f_calc_Req(R0, bubtherm, masstrans, ee, vmaterial);
    %     % initial velocity
    %     Uzero = -(1-pzero)/(Cstar); 
    % else
    %     % Plesset & Prosperetti, ARFM 1977, p166
    %     Uzero = 0;
    % end
    % 
    % % Plesset & Prosperetti, ARFM 1977, p166
    % if  (wave_type == -2)
    %     % initial velocity
    %     Uzero = -(ee/P8)/(Cstar); 
    % end
    
    % out parameters

    % equation settings
    eqns_opts = [radial bubtherm medtherm stress eps3 vapor masstrans];
    % solver options
    solve_opts = [method spectral divisions Nv Nt Mt Lv Lt];
    % dimensionless initial conditions
    init_opts = [Rzero Uzero P0_star P8 T8 Pv_star Req_zero alphax];
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