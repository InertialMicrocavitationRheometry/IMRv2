% file m_imr_spectral.m
% brief contains module m_imr_spectral

% brief This module features a Chebyshev spectral collocation solver of the
% PDEs involving thermal transport and viscoelasticity to solve
% Rayleigh-Plesset equations
function varargout = f_imr_spectral(varargin)
    
    % problem initialization
    [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
        acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] ...
        = f_call_params(varargin{:});
    
    % defaults
    
    % non-Newtonian viscosity parameters
    fnu             = 0;
    intfnu          = 0;
    dintfnu         = 0;
    ddintfnu        = 0;
    
    % equations settings
    radial          = eqns_opts(1);
    bubtherm        = eqns_opts(2);
    medtherm        = eqns_opts(3);
    stress          = eqns_opts(4);
    eps3            = eqns_opts(5);
    masstrans       = eqns_opts(6);
    % if (stress == 4)
    %     ptt = 1;
    % else
    %     ptt = 0;
    % end
    
    % solver options
    method          = solve_opts(1);
    spectral        = solve_opts(2);
    divisions       = solve_opts(3);
    Nv              = solve_opts(4);
    Nt              = solve_opts(5);
    Mt              = solve_opts(6);
    Lv              = solve_opts(7);
    Lt              = solve_opts(8);
    
    % dimensionless initial conditions
    Rzero           = init_opts(1);
    Uzero           = init_opts(2);
    Pb_star         = init_opts(3);
    % P8              = init_opts(4);
    T8              = init_opts(5);
    Pv_star         = init_opts(6);
    Req             = init_opts(7);
    alphax          = init_opts(8);
    
    % dimensionaless initial stress
    Szero           = init_stress;
    
    % time span options
    tspan = tspan_opts;
    
    % output options
    dimensionalout  = out_opts(1);
    progdisplay     = out_opts(2);
    tref            = out_opts(3);
    Rref            = out_opts(4);
    Uref            = out_opts(5);
    
    % physical parameters
    
    % acoustic parameters
    Cstar           = acos_opts(1);
    GAMa            = acos_opts(2);
    kappa           = acos_opts(3);
    nstate          = acos_opts(4);
    hugoniot_s      = acos_opts(5);
    
    % dimensionless waveform parameters
    om              = wave_opts(1);
    ee              = wave_opts(2);
    tw              = wave_opts(3);
    dt              = wave_opts(4);
    mn              = wave_opts(5);
    wave_type       = wave_opts(6);
    if wave_type < 0
        wave_poly   = wave_opts(7);
        wave_dpoly  = wave_opts(8);
    else
        wave_poly   = [];
        wave_dpoly  = [];
    end
    
    pvarargin = [om,ee,tw,dt,mn,wave_type,wave_poly,wave_dpoly];
    
    % dimensionless viscoelastic
    We              = sigma_opts(1);
    Re8             = sigma_opts(2);
    DRe             = sigma_opts(3);
    v_a             = sigma_opts(4);
    v_nc            = sigma_opts(5);
    Ca              = sigma_opts(6);
    LAM             = sigma_opts(7);
    De              = sigma_opts(8);
    JdotA           = sigma_opts(9);
    nu_model        = sigma_opts(10);
    v_lambda_star   = sigma_opts(11);
    zeNO            = sigma_opts(12);
    iDRe            = sigma_opts(13);
    iWe             = 1/We;
    
    % dimensionless thermal
    Foh             = thermal_opts(1);
    Br              = thermal_opts(2);
    alpha_g         = thermal_opts(3);
    beta_g          = thermal_opts(4);
    % alpha_v         = thermal_opts(5);
    % beta_v          = thermal_opts(6);
    chi             = thermal_opts(7);
    iota            = thermal_opts(8);
    
    % dimensionaless mass transfer
    Fom             = mass_opts(1);
    iota = Fom*0 + iota;
    % C0              = mass_opts(2);
    % Rv_star         = mass_opts(3);
    % Ra_star         = mass_opts(4);
    % L_heat_star     = mass_opts(5);
    % mv0             = mass_opts(6);
    % ma0             = mass_opts(7);
    
    % pre_process code
    
    % collocation point construction
    y = cos(pi*(0:Nt)'/(2*Nt));
    xi = cos(pi*(0:Mt)'/Mt);
    % ze = cos(pi*(1:Nv)'/Nv);
    % collocation matrix construction
    [gA,gAI,~,~,gAPd,gAPdd] = f_dcdmtxe(Nt);
    [mA,~,~,~,mAPd,mAPdd] = f_dcdmtx(Mt);
    % [gC,gCI,~,~,~,~] = f_dcdmtxe(Nt);
    Q = [gA(2:end,:) zeros(Nt,Mt+1);
    zeros(Mt,Nt+1) mA(2:end,:);
    2*(0:Nt).^2 iota*(0:Mt).^2];
    Q = sparse(Q);
    % [sCA,sCI,sCAd,~,~,~] = f_dcdmtx(Nv);
    % sCA = sCA(2:end,2:end) - 1;
    % sCI = sCI(2:end,2:end);
    % sCAd = sCAd(2:end,2:end);
    
    % precomputations
    % LDR = LAM*De/Re8;
    sam = 1 + GAMa;
    no = (nstate-1)/nstate;
    nog = (nstate-1)/2;
    kapover = (kappa-1)/kappa;
    yT = 2*Lt./(1+xi) - Lt + 1;
    yT2 = yT.^2;
    yT3 = yT.^3;
    iyT3 = yT.^-3;
    iyT4 = yT.^-4;
    iyT6 = yT.^-6;
    % yV = 2*Lv./(1-ze) - Lv + 1;
    nn = ((-1).^(0:Nt).*(0:Nt).^2)';
    nn = sparse(nn);
    
    % precomputations for viscous dissipation
    zT = 1 - 2./(1 + (yT - 1)/Lv);
    ZZT = cos(acos(zT)*(1:Nv)) - 1;
    cdd = f_preStressInt(Lv,Nv);
    
    % index management
    if spectral == 0
        Nv = 1;
    end
    if bubtherm == 0
        Nt = -1;
        Mt = -1;
        qdot = [];
    end
    if medtherm == 0
        Mt = -1;
    end
    % TODO add masstransfer back
    % if masstrans == 0
    %     Nm = -1;
    % end
    ia = 4:(4+Nt);
    ib = (5+Nt):(5+Nt+Mt);
    ivisco1 = (6+Nt+Mt):(5+Nt+Mt+Nv);
    ivisco2 = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);
    % ie = (6+Nt+Mt+2*Nv):(5+Nt+Mt+2*Nv+Nm);
    
    % initial condition assembly
    
    % radius, velocity, pressure
    
    % auxiliary temperature, boundary temperature, medium temperature,
    Tau0 = zeros(Nt+1,1);
    Tm0 = ones(Mt ~= -1);
    Tm1 = zeros(Mt,1);
    
    % TODO ADD THE MASS Transfer structure here
    
    % initial condition vector
    init = [Rzero;
    Uzero;
    Pb_star;
    Tau0;
    Tm0;
    Tm1;
    Szero];
    
    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, 'spectral');
    bubble = @SVBDODE;
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan);
    
    % post processing
    
    % extract result
    R = X(:,1);
    Rdot = X(:,2);
    P = X(:,3);
    % extracting the Chebyshev coefficients
    a = X(:,ia)';
    b = X(:,ib)';
    if bubtherm
        T = (alpha_g-1+sqrt(1+2*alpha_g*gA*a))/alpha_g;
        if medtherm
            Tm = mA*b;
        end
    else
        T = R.^(-3*kappa);
    end
    % Z1 = X(:,ivisco1);
    % Z2 = X(:,ivisco2);
    % if masstrans
    %     C = gC*e;
    % end
    
    pA = zeros(size(t));
    for n = 1:length(t)
        pA(n) = f_pinfinity(t(n),pvarargin);
    end
    
    % transform to real space
    if spectral == 1
        % trr = sCA*c;
        % t00 = sCA*d;
    else
        %trr = c;
        %t00 = d;
    end
    
    % transform variables back into their dimensional form
    if dimensionalout == 1
        t = t*tref;
        R = R*Rref;
        Rdot = Rdot*Uref;
        P = P*P8;
        T = T*T8;
        %pA = pA*p0;
        % c = c*p0;
        % d = d*p0;
        % e = e*C0;
        % C = C*C0;
        Rddot = Rddot*uc/tc;
        if spectral == 1
            % trr = trr*p0;
            % t00 = t00*p0;
        end
        if bubtherm == 1
            if medtherm == 1
                Tm = Tm*T8;
            end
        end
    end
    
    % outputs
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = Rdot;
    varargout{4} = P;
    if bubtherm == 1
        varargout{5} = T;
    else
        varargout{5} = [];
    end
    if medtherm == 1
        varargout{6} = Tm;
    else
        varargout{6} = ((T8 - 1)*dimensionalout + 1)*ones(size(t,1),1);
    end
    varargout{7} = [];
    
    % solver function
    function dXdt = SVBDODE(t,X)
        
        % showing output
        if progdisplay
            disp(t/tfin);
        end
        
        % extract standard inputs
        
        % bubble wall radius
        R = X(1);
        % bubble wall velocity
        Rdot = X(2);
        % internal bubble pressure
        P = X(3);
        % thermal vector in and outside the bubble
        qdot = [];
        
        % updating the viscous forces/Reynolds number
        if nu_model ~= 0
            [fnu,intfnu,dintfnu,ddintfnu] = f_viscosity(nu_model,Rdot, ...
                R,v_a,v_nc,v_lambda_star);
        end
        
        % bubble equations of motion
        if bubtherm
            % extract auxiliary temperature
            theta = gA*X(ia);
            % auxiliary temperature derivatives
            
            % first order derivative
            dtheta = gAPd*theta;
            % second order derivative
            ddtheta = gAPdd*theta;
            % temperature and thermal diffusivity fields
            T = (alpha_g - 1 + sqrt(1+2*alpha_g*theta))/alpha_g;
            D = kapover*(alpha_g*T.^2 + beta_g*T)/P;
            
            % bubble pressure equation
            Pdot = 3/R*((kappa-1)*chi/R*dtheta(1) - kappa*P*Rdot);
            
            % bubble wall velocity equation
            U_vel = chi/R^2*(2*D./y - kapover/P*(dtheta-dtheta(1)*y));
            
            % auxiliary temperature derivative
            thetadot = Pdot*D + U_vel.*dtheta + chi*D/R^2.*ddtheta;
            
            thetadot(end) = Pdot*D(end) - chi/R^2*(8*D(end)*sum(nn.*X(ia)) ...
                + kapover/P*dtheta(end)^2) + chi*D(end)/R^2.*ddtheta(end);
            
            % surroundings equations of motion
            
            % surrounding temperature
            if medtherm
                % extract medium temperature
                Tm = mA*X(ib);
                
                % material temperature equations
                advection = (1+xi).^2/(Lt*R).*...
                    (Foh/R*((1+xi)/(2*Lt) - 1./yT) + ...
                    Rdot/2*(1./yT2 - yT)).*(mAPd*Tm);
                diffusion = Foh/4*(1+xi).^4/(Lt^2*R^2).*(mAPdd*Tm);
                % stress dissipation
                [taugradu] = f_stress_dissipation(stress,spectral,Req,R,Rdot, ...
                    Ca,Br,Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1, ...
                    ivisco2,fnu,DRe);
                % surrounding temperature evolution equation
                Tmdot = advection + diffusion + taugradu;
                
                % enforce boundary condition and solve A*x=b problem
                Tmdot(end) = 0;
                qdot = [ones(1,Nt+1) -(alpha_g*(T(1)-1)+1)*ones(1,Mt+1); Q]...
                    \[0;
                thetadot(2:end);
                Tmdot(2:end);
                0];
            else
                % cold-liquid approximation
                % solve auxiliary temperature with boundary condition
                qdot = gAI*[0;
                thetadot(2:end)];
            end
        else
            % polytropic gas + vapor
            P = (Pb_star-Pv_star)*(1/R)^(3*kappa) + Pv_star;
            Pdot = -3*kappa*P*Rdot/R;
        end
        
        % pressure waveform
        [Pf8,Pf8dot] = f_pinfinity(t,pvarargin);
        
        % stress equation evolution
        [S,Sdot,Z1dot,Z2dot] = f_stress_calc(stress,X,Req,R,Ca,De,Re8, ...
            Rdot,alphax,ivisco1,ivisco2,LAM,zeNO,cdd,intfnu,dintfnu,iDRe);
        
        % bubble wall evolution / acceleration
        [Rddot] = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, S, ...
            Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
            ddintfnu, iDRe);
        
        % output assembly
        dXdt = [Rdot;
        Rddot;
        Pdot;
        qdot;
        Z1dot;
        Z2dot];
        
    end
    % end of solver
    
    % function Cw= CW(Tw,P)
    %   % Calculates the concentration at the bubble wall
    %   %Function of P and temp at the wall
    %   thetha = Rv_star/Ra_star*(P./(f_pvsat(Tw*T8)/P8) -1);
    %   Cw = 1./(1+thetha);
    % end
    
end

% precomputation functions

% discrete Chebyshev derivative matrices, FCD
function [A,B,D,E,C,F] = f_dcdmtx(N)
    
    jj = 0:N;
    jj2 = jj.^2;
    theta = pi*jj'/N;
    
    % matrix for a -> p
    A = cos(theta*jj);
    
    % matrix for p -> a
    B = 2*cos(theta*jj)'/N;
    B(:,[1 N+1]) = 0.5*B(:,[1 N+1]);
    B([1 N+1],:) = 0.5*B([1 N+1],:);
    
    theta = theta(2:N);
    
    % matrix for a -> dp/dx
    D = zeros(N+1);
    D(1,:) = jj2;
    D(2:N,:) = csc(theta)*jj.*sin(theta*jj);
    D(N+1,:) = jj2.*((-1).^(jj+1));
    
    % matrix for a -> d2p/dx2
    E = zeros(N+1);
    E(1,:) = jj2.*(jj2-1)/3;
    E(2:N,:) = (cos(theta)./sin(theta).^3)*jj.*sin(theta*jj) ...
        - (csc(theta).^2)*((0:N).^2).*cos(theta*jj);
    E(N+1,:) = jj2.*(jj2-1).*(-1).^jj/3;
    
    % matrix for p -> dp/dx
    C = D*B;
    
    % matrix for p -> d2p/dx2
    F = E*B;
    
end

% even discrete Chebyshev derivative matrices, FCD
function [A,B,D,E,C,F] = f_dcdmtxe(N)
    
    jj = 2*(0:N);
    jj2 = jj.^2;
    theta = pi*jj'/(4*N);
    
    % matrix for a -> p
    A = cos(theta*jj);
    
    % matrix for p -> a
    B = 2*cos(theta*jj)'/N;
    B(:,[1 N+1]) = 0.5*B(:,[1 N+1]);
    B([1 N+1],:) = 0.5*B([1 N+1],:);
    
    theta = theta(2:N+1);
    
    % matrix for a -> dp/dx
    D = zeros(N+1);
    D(1,:) = jj2;
    D(2:N+1,:) = csc(theta)*jj.*sin(theta*jj);
    
    % matrix for a -> d2p/dx2
    E = zeros(N+1);
    E(1,:) = jj2.*(jj2-1)/3;
    E(2:N+1,:) = (cos(theta)./sin(theta).^3)*jj.*sin(theta*jj) ...
        - (csc(theta).^2)*(jj2).*cos(theta*jj);
    
    % matrix for p -> dp/dx
    C = D*B;
    
    % matrix for p -> d2p/dx2
    F = E*B;
    
end

% prestress integrals
function cdd = f_preStressInt(L,N)
    
    % integral precomputations
    Lstr = ['L' num2str(L,18)];
    Lstr = strrep(Lstr,'.','p');
    if exist('d_stress_store.mat','file') ~= 0
        load('d_stress_store.mat','store');
        if isfield(store,Lstr) == 1
            if size(store.(Lstr),2) >= N, Nstart = 0;
            else
                Nstart = size(store.(Lstr),2) + 1;
                disp('Past integral precomputation not found in full, catching up ...');
            end
        else
            Nstart = 1;
            disp('Past integral precomputation not found, starting anew ...');
        end
    else
        Nstart = 1;
        store = struct;
        disp('Past integral precomputation not found, starting anew ...');
    end
    if Nstart ~= 0 % begin extended precomputation
        
        store.(Lstr)(Nstart:N) = f_StressInt(L,N,Nstart);
        save('stress_store.mat','store');
        disp('Precomputation completed.');
        
    end
    cdd = store.(Lstr)(1:N)';
    
end

%
function cdd = f_StressInt(L,N,varargin)
    
    if nargin == 2
        k = 1;
    else
        k = varargin{1};
    end
    syms x;
    cdd = zeros(N-k+1,1);
    
    for n = k:N
        cdd(n-k+1) = subs(2*L*int((cos(n*acos(x))-1)/((L*(2/(1-x)-1)+1)*(1-x)^2),-1,1));
    end
    
end
