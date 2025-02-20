% file m_imrv2_spectral.m
% brief contains module m_imrv2_spectral

% brief This module features a Chebyshev spectral collocation solver of the
% PDEs involving thermal transport and viscoelasticity to solve
% Rayleigh-Plesset equations
function varargout =  m_imrv2_spectral(varargin)
    
    % problem Initialization
    [eqns_opts, solve_opts, init_opts, tspan_opts, out_opts, acos_opts,...
        wave_opts, sigma_opts, thermal_opts, mass_opts]...
    = f_call_params(varargin{:});
    
    % equations settings
    radial          = eqns_opts(1);
    bubtherm        = eqns_opts(2);
    medtherm        = eqns_opts(3);
    stress          = eqns_opts(4);
    eps3            = eqns_opts(5);
    vapor           = eqns_opts(6);
    masstrans       = eqns_opts(7);
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
    p0star          = init_opts(3);
    P8              = init_opts(4);
    T8              = init_opts(5);
    Pv_star         = init_opts(6);
    Req             = init_opts(7);
    alphax          = init_opts(8);
    
    % time span options
    tspan = tspan_opts;
    tfin = tspan(end);
    % output options
    dimensionalout  = out_opts(1);
    progdisplay     = out_opts(2);
    
    % physical parameters
    
    % acoustic parameters
    Cstar           = acos_opts(1);
    GAMa            = acos_opts(2);
    kappa           = acos_opts(3);
    nstate          = acos_opts(4);
    % dimensionless waveform parameters
    om              = wave_opts(1);
    ee              = wave_opts(2);
    tw              = wave_opts(3);
    dt              = wave_opts(4);
    mn              = wave_opts(5);
    wave_type       = wave_opts(6);
    
    pvarargin = [om,ee,tw,dt,mn,wave_type];
    
    % dimensionless viscoelastic
    We              = sigma_opts(1);
    Re8             = sigma_opts(2);
    % DRe             = sigma_opts(3);
    % v_a             = sigma_opts(4);
    % v_nc            = sigma_opts(5);
    Ca              = sigma_opts(6);
    LAM             = sigma_opts(7);
    De              = sigma_opts(8);
    JdotA           = sigma_opts(9);
    % vmaterial       = sigma_opts(10);
    % v_lambda_star   = sigma_opts(11);
    zeNO            = sigma_opts(12);
    iWe             = 1/We;
    if Ca==-1
        Ca=Inf;
    end
    % dimensionless thermal
    Foh             = thermal_opts(1);
    Br              = thermal_opts(2);
    alpha           = thermal_opts(3);
    beta            = thermal_opts(4);
    chi             = thermal_opts(5);
    iota            = thermal_opts(6);
    % dimensionaless mass transfer
    % Fom             = mass_opts(1);
    C0              = mass_opts(2);
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
    [gA,gAI,~,~,gAPd,gAPdd] = dcdmtxe(Nt);
    [mA,~,~,~,mAPd,mAPdd] = dcdmtx(Mt);
    % [gC,gCI,~,~,~,~] = dcdmtxe(Nt);
    Q = [gA(2:end,:) zeros(Nt,Mt+1);
    zeros(Mt,Nt+1) mA(2:end,:);
    2*(0:Nt).^2 iota*(0:Mt).^2];
    Q = sparse(Q);
    % [sCA,sCI,sCAd,~,~,~] = dcdmtx(Nv);
    % sCA = sCA(2:end,2:end) - 1;
    % sCI = sCI(2:end,2:end);
    % sCAd = sCAd(2:end,2:end);
    
    % precomputations
    % LDR = LAM*De/Re8;
    sam = 1 - Pv_star + GAMa;
    no = (nstate-1)/nstate;
    kapover = (kappa-1)/kappa;
    yT = 2*Lt./(1+xi) - Lt + 1;
    % yV = 2*Lv./(1-ze) - Lv + 1;
    nn = ((-1).^(0:Nt).*(0:Nt).^2)';
    nn = sparse(nn);
    Udot = 0*C0; %TODO Change this!
    
    % precomputations for viscous dissipation
    zT = 1 - 2./(1 + (yT - 1)/Lv);
    ZZT = cos(acos(zT)*(1:Nv)) - 1;
    cdd = preStressInt(Lv,Nv);
    
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
    if masstrans == 0
        Nm = -1;
    end
    ia = 4:(4+Nt);
    ib = (5+Nt):(5+Nt+Mt);
    ic = (6+Nt+Mt):(5+Nt+Mt+Nv);
    id = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);
    ie = (6+Nt+Mt+2*Nv):(5+Nt+Mt+2*Nv+Nm);
    
    % initial condition assembly
    
    % radius, velocity, pressure, auxiliary temperature spectrum,
    % medium temperature spectrum, stress spectrum,
    % initial stress integral
    Tau0 = zeros(Nt+1,1);
    Tm0 = ones(Mt ~= -1);
    Tm1 = zeros(Mt,1);
    % TODO ADD THE MASS Transfer structure here
    
    Sp = zeros(2*(Nv - 1)*(spectral == 1) + 2,1);
    init = [Rzero;
    Uzero;
    p0star;
    Tau0;
    Tm0;
    Tm1;
    Sp;
    0];
    
    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        eps3, Re8, De, Ca, LAM, 'spectral');
    stepcount = 0;
    bubble = @SVBDODE;
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);
    
    % post processing
    
    % extract result
    R = X(:,1);
    U = X(:,2);
    p = X(:,3);
    % extracting the Chebyshev coefficients
    Z1 = X(:,ic);
    Z2 = X(:,id);
    a = X(:,ia)';
    b = X(:,ib)';
    c = X(:,ic)';
    d = X(:,id)';
    % e = X(:,ie)';
    I = X(:,ie+1);
    if bubtherm
        T = (alpha-1+sqrt(1+2*alpha*gA*a))/alpha;
        if medtherm
            TL = mA*b;
        end
    else
        T = R.^(-3*kappa);
    end
    % if masstrans
    %     C = gC*e;
    % end
    
    pA = zeros(size(t));
    for n = 1:length(t)
        pA(n) = f_pinfinity(t(n),pvarargin);
    end
    
    % transform to real space
    if spectral == 1
        trr = sCA*c;
        t00 = sCA*d;
    else
        trr = c;
        t00 = d;
    end
    % dimensionalization
    if dimensionalout == 1
        % re-dimensionalize problem
        t = t*tc;
        R = R*Rref;
        U = U*uc;
        p = p*p0;
        T = T*T8;
        %pA = pA*p0;
        I = I*p0;
        % c = c*p0;
        % d = d*p0;
        % e = e*C0;
        % C = C*C0;
        Udot = Udot*uc/tc;
        if spectral == 1
            trr = trr*p0;
            t00 = t00*p0;
        end
        if bubtherm == 1
            if medtherm == 1
                TL = TL*T8;
            end
        end
    end
    
    % outputs
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = U;
    varargout{4} = p;
    varargout{5} = T;
    if bubtherm == 1 && medtherm == 1
        varargout{6} = TL;
    else
        varargout{6} = ((T8 - 1)*dimensionalout + 1)*ones(size(t,1),1);
    end
    if masstrans == 1
        varargout{7} = C;
    end
    varargout{8} = trr;
    varargout{9} = t00;
    varargout{10} = I;
    if stress == 2
        varargout{11} = Z1; % Z1
        varargout{12} = Z2; % Z2
        varargout{13} = (Z1)./R.^3 - 4*LAM/Re8.*U./R; % J
    elseif stress == 3
        varargout{11} = Z1; % Z1
        varargout{12} = Z2; % Z2
        varargout{13} = (Z1 + Z2)./R.^3 - 4*LAM./Re8.*U./R; % J
    end
    
    % solver function
    function dXdt = SVBDODE(t,X)
        stepcount = stepcount + 1;
        if progdisplay == 1, disp(t/tfin); end
            
            % extract standard inputs
            R = X(1);
            U = X(2);
            p = X(3);
            qdot = [];
            
            % updating the viscous forces/Reynolds number
            % [fnu,intfnu,dintfnu,ddintfnu] = ...
                % [fnu,~,~,~] = ...
            % f_nonNewtonian_integrals(vmaterial,U,R,v_a,v_nc,v_lambda_star);
            
            % non-condensible gas pressure and temperature
            if bubtherm
                % extract auxiliary temperature
                SI = gA*X(ia);
                % auxiliary temperature derivatives
                
                % first order derivative
                dSI = gAPd*SI;
                % second order derivative
                ddSI = gAPdd*SI;
                % temperature and thermal diffusivity fields
                T = (alpha - 1 + sqrt(1+2*alpha*SI))/alpha;
                D = kapover*(alpha*T.^2 + beta*T)/p;
                
                % vapor pressure
                pVap = vapor*f_pvsat(T(end)*T8)/P8;
                
                % bubble pressure
                pdot = 3/R*((kappa-1)*chi/R*dSI(1) - kappa*p*U);
                
                % auxiliary temperature derivative
                SIdot = pdot*D + chi/R^2*(2*D./y - kapover/p*(dSI - dSI(1)*y)).*dSI ...
                    + chi*D/R^2.*ddSI;
                
                SIdot(end) = pdot*D(end) - chi/R^2*(8*D(end)*sum(nn.*X(ia)) ...
                    + kapover/p*dSI(end)^2) + chi*D(end)/R^2.*ddSI(end);
                
                if medtherm % warm-liquid
                    % extract medium temperature
                    TL = mA*X(ib);
                    % new derivative of medium temperature
                    first_term = (1+xi).^2/(Lt*R).*...
                        (Foh/R*((1+xi)/(2*Lt) - 1./yT) + ...
                    U/2*(1./yT.^2 - yT)).*(mAPd*TL);
                    second_term = Foh/4*(1+xi).^4/(Lt^2*R^2).*(mAPdd*TL);
                    TLdot =  first_term + second_term;
                    % include viscous heating
                    if spectral == 1
                        TLdot = TLdot - 2*Br*U./(R*yT.^3).*(ZZT*(X(ic) - X(id)));
                    else %if stress == 1
                        TLdot = TLdot + 4*Br./yT.^6*(U/R*(1-1/R^3)/Ca + 3/Re8*(U/R)^2);
                    end
                    % enforce boundary condition and solve
                    TLdot(end) = 0;
                    qdot = [ones(1,Nt+1) -(alpha*(T(1)-1)+1)*ones(1,Mt+1); Q]...
                        \[0; SIdot(2:end); TLdot(2:end); 0];
                else % cold-liquid approximation
                    % solve auxiliary temperature with boundary condition
                    qdot = gAI*[0;
                    SIdot(2:end)];
                end
            else
                % polytropic approximation
                pVap = vapor*Pv_star;                
                pdot = -3*kappa*U/R*p;
            end
            
            % stress equation
            [J,JdotX,Z1dot,Z2dot] = ...
                f_stress_calc(stress,X,Req,R,Ca,De,Re8,U,alphax,ic,id,LAM,zeNO,cdd);
            
            % pressure waveform
            [pf8,pf8dot] = f_pinfinity(t,pvarargin);
            
            % bubble wall acceleration
            [Udot] = f_radial_eq(radial, p, pdot, pVap, pf8, pf8dot, iWe, R, U, ...
                J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA );
            
            % stress integral rate
            Jdot = JdotX - JdotA*Udot/R;
            
            % output assembly
            dXdt = [U;
            Udot;
            pdot;
            qdot;
            Z1dot;
            Z2dot;
            Jdot];
            
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
    function [A,B,D,E,C,F] = dcdmtx(N)
        %FCD	Discrete Chebyshev derivative matrices
        
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
    
    function [A,B,D,E,C,F] = dcdmtxe(N)
        %FCD	Even discrete Chebyshev derivative matrices
        
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
    
    
    function cdd = preStressInt(L,N)
        
        % integral precomputations
        Lstr = ['L' num2str(L,18)];
        Lstr = strrep(Lstr,'.','p');
        if exist('StressIntStore.mat','file') ~= 0
            load('StressIntStore.mat','store');
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
            
            store.(Lstr)(Nstart:N) = StressInt(L,N,Nstart);
            save('StressIntStore.mat','store');
            disp('Precomputation completed.');
            
        end
        cdd = store.(Lstr)(1:N)';
        
    end
    
    function cdd = StressInt(L,N,varargin)
        
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
