% file m_imrv2_finitediff.m
% brief contains module m_imrv2_finitediff

% brief This module features a fourth- and sixth-order accurate finite
% difference solver of the PDEs involving thermal transport and
% viscoelasticity to solve Rayleigh-Plesset equations
function varargout =  m_imrv2_finitediff(varargin)
    
    % problem initialization
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
    Fom             = mass_opts(1);
    C0              = mass_opts(2);
    Rv_star         = mass_opts(3);
    Ra_star         = mass_opts(4);
    L_heat_star     = mass_opts(5);
    % mv0             = mass_opts(6);
    % ma0             = mass_opts(7);
    
    % pre_process
    
    % creates finite difference matrices
    D_Matrix_T_C = f_finite_diff_mat(Nt,1,0);
    DD_Matrix_T_C = f_finite_diff_mat(Nt,2,0);
    D_Matrix_Tm = f_finite_diff_mat(Mt,1,1);
    DD_Matrix_Tm = f_finite_diff_mat(Mt,2,1);
    
    % create spatial nodes
    
    % inside the bubble
    N = Nt-1;
    deltaY = 1/N;
    i = 1:1:N+1;
    y = ((i-1)*deltaY)';
    % outside the bubble
    Nm = Mt-1;
    deltaYm = -2/Nm;
    j = 1:1:Nm+1;
    
    % precomputations
    %LDR = LAM*De/Re8;
    sam = 1 - Pv_star + GAMa;
    no = (nstate-1)/nstate;
    kapover = (kappa-1)/kappa;
    xi = (1+(j-1)*deltaYm)';
    yT = ((2./(xi+1)-1)*Lt+1);
    
    % index management
    if spectral == 0
        Nv = 0;
    end
    if bubtherm == 0
        Nt = 0;
        Mt = 0;
    end
    if medtherm == 0
        Mt = 0;
    end
    ic = (4+Nt+Mt):(4+Nt+Mt+Nv);
    id = (5+Nt+Mt+Nv):(4+Nt+Mt+2*Nv);
    
    % precomputations for viscous dissipation
    % zT = 1 - 2./(1 + (yT - 1)/Lv);
    cdd = preStressInt(Lv,Nv);
    
    % initial condition assembly
    
    % radius, velocity, pressure, bubble temperature, medium temperature,
    % vapor concentration
    T = zeros(-1,1);
    % bubble temperature initial condition
    if bubtherm
        Tau0 = zeros(Nt,1);
    else
        Tau0 = zeros(-1,1);
    end
    % medium initial temperature
    if medtherm
        Tm0 = ones(Mt,1);
    else
        Tm0 = zeros(-1,1);
    end
    % mass transfer initial condition
    if masstrans
        C0 = C0*ones(Nt,1);
    else
        C0 = zeros(-1,1);
    end
    % stress spectra
    if stress < 3
        istress = 0;
    elseif stress == 3 || stress == 4
        istress = 1;
    elseif stress == 5
        istress = 2;
    end
    Sp = zeros(2*(Nv - 1)*(spectral == 1) + istress,1);
    
    % initial condition vector
    init = [Rzero;
    Uzero;
    p0star;
    Tau0;
    Tm0;
    C0;
    Sp];
    
    % thermal auxiliary variable for boundary conditions
    tau_del = [];
    TL = [];
    
    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral, ...
        eps3, Re8, De, Ca, LAM, 'finite difference');
    bubble = @SVBDODE;
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);
    
    % post processing
    
    % extract result
    R = X(:,1);
    U = X(:,2);
    p = X(:,3);
    if bubtherm
        Tau = X(:,4:(Nt+3));
        T = (alpha - 1 + sqrt(1+2*Tau*alpha)) / alpha; % Temp in bubble
        if medtherm
            Tm = X(:,(Nt+4):(2*Nt+3));
        end
    end
    if masstrans
        C = X(:,(2*Nt+4):end);
    end
    
    % transform variables back into their dimensional form
    if (dimensionalout == 1)
        t = t*tc;
        R = R*Rref;
        U = U*uc;
        p = p*p0;
        if bubtherm == 1
            T = T*T8;
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
    
    % TODO Add the stress output if possible, similar to the spectral code
    
    % solver function
    function dXdt = SVBDODE(t,X)
        if progdisplay == 1
            disp(t/tfin);
        end
        
        % extract standard inputs
        R = X(1);
        U = X(2);
        p = X(3);
        
        % solve for boundary condition at the wall
        if (medtherm)
            Tau = X(4:(Nt+3));
            Tm = X((Mt+4):(2*Mt+3));
            if t/tfin > 0.001
                %Might need to tune 0.001 for convergence:
                guess= -.001+tau_del(end);
                prelim  = fzero(@Boundary,guess);
            else
                guess = -.0001;
                prelim  = fzero(@Boundary,guess);
            end
        else
            prelim = 0 ;
        end
        
        % sets value at boundary conditions
        if bubtherm
            % extract the temperature
            Tau = X(4:(Nt+3));
            % using the previous bubble wall temperature as an initial guess for speed up
            tau_del=[tau_del prelim];
            % setting the temperature at the bubble wall
            Tau(end) = prelim;
            T = TW(Tau);
            
            % precalculation for speed up
            K_star = alpha*T+beta;
        end
        if medtherm
            Tm(1) = T(end);
        end
        
        % updating the viscous forces/Reynolds number
        % [fnu,intfnu,dintfnu,ddintfnu] = ...
            % [fnu,~,~,~] = ...
        % f_nonNewtonian_integrals(vmaterial,U,R,v_a,v_nc,v_lambda_star);
        
        Taudot = zeros(-1,1);
        Tmdot = zeros(-1,1);
        Cdot = zeros(-1,1);
        % heat and mass transfer PDE residual calculations
        if bubtherm && masstrans
            C = X((2*Nt+4):end);
            C(end) =  CW(T(end),P);
            Rmix = C*Rv_star + (1-C)*Ra_star;
            % concentration field inside the bubble
            DC  = D_Matrix_T_C*C;
            DDC = DD_Matrix_T_C*C;
            % temp. field inside the bubble
            DTau  = D_Matrix_T_C*Tau;
            DDTau = DD_Matrix_T_C*Tau;
            
            % internal pressure equation
            pVap = vapor*(f_pvsat(T(1)*T8)/P8);
            pdot = 3/R*(chi*(kappa-1)*DTau(end)/R-kappa*p*U+...
                + kappa*p*Fom*Rv_star*DC(end)...
            /( T(end)*R*Rmix(end)*(1-C(end)) ) );
            
            % temperature inside the bubble
            U_vel = (chi/R*(kappa-1)*DTau-y*R*pdot/3)/(kappa*p);
            first_term = (DDTau*chi./R^2+pdot)*(K_star*T./p*kapover);
            second_term = -DTau.*((1./R).*(U_vel-y.*U));
            
            Taudot= first_term+second_term;
            Taudot(end) = 0;
            
            % vapor concentration equation
            U_mix = U_vel + Fom/R*((Rv_star - Ra_star)/Rmix)*DC  ;
            one = DDC;
            two = DC*(DTau/(K_star*T)+((Rv_star - Ra_star)/Rmix)*DC );
            three =  (U_mix-U*yk)/R*DC;
            
            % concentration evolution
            Cdot = Fom/R^2*(one - two) - three;
            Cdot(end) = 0;
            
        elseif bubtherm
            % temp. field inside the bubble
            DTau  = D_Matrix_T_C*Tau;
            DDTau = DD_Matrix_T_C*Tau;
            
            % internal pressure equation
            pVap = vapor*(f_pvsat(T(1)*T8)/P8);
            pdot = 3/R*(chi*(kappa-1)*DTau(end)/R-kappa*p*U);
            
            % temperature inside the bubble
            
            % first_term = (DDTau.*chi./R^2+pdot).*( K_star.*T/p*(kappa-1)/kappa);
            first_term = (DDTau.*chi./R^2+pdot).*(kapover*K_star.*T./p);
            U_vel = (chi/R*(kappa-1)*DTau-y*R*pdot/3)/(kappa*p);
            second_term = -DTau.*((1/R).*(U_vel-y*U));
            
            Taudot= first_term+second_term;
            Taudot(end) = 0;
        else
            % polytropic gas
            pVap = vapor*Pv_star;
            pdot= -3*kappa*U/R*p;
        end
        
        if medtherm
            % temp. field outside the bubble
            DTm = D_Matrix_Tm*Tm;
            DDTm = DD_Matrix_Tm*Tm;
            % warm liquid
            first_term = (1+xi).^2./(Lt*R).*...
                (U./yT.^2.*(1-yT.^3)/2+Foh/R.*((xi+1)/(2*Lt)-1./yT)).*DTm;
            % second_term = foh/R^2.*(xk+1).^4/L^2.*DDTm/4;
            second_term = Foh/R^2.*(xi+1).^4/Lt^2.*DDTm/4;
            third_term =  4*Br./yT.^6.*(3/Re8.*(U/R)^2);
            %third_term =  4*Br./yT.^6.*(3/(Re8+DRe*fnu).*(U/R)^2);
            Tmdot = first_term+second_term+third_term;
            % Sets boundary condition on temp
            Tmdot(end) = 0;
            % Previously calculated;
            Tmdot(1) = 0;
        end
        
        % stress equation
        [J,JdotX,Z1dot,Z2dot] = ...
            f_stress_calc(stress,X,Req,R,Ca,De,Re8,U,alphax,ic,id,LAM,zeNO,cdd);
        
        % pressure waveform
        [pf8,pf8dot] = f_pinfinity(t,pvarargin);
        
        % bubble wall acceleration
        [Udot] = f_radial_eq(radial, p, pdot, pVap, pf8, pf8dot, iWe, R, U, ...
            J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA );
        
        % output assembly
        dXdt = [U;
        Udot;
        pdot;
        Taudot;
        Tmdot;
        Cdot;
        Z1dot;
        Z2dot];
        
    end
    % end of solver
    
    % functions called by solver
    
    function Tw= TW(Tauw)
        %calculates the temperature at the bubble wall as a function of \tau
        Tw = (alpha - 1 + sqrt(1+2*Tauw*alpha)) / alpha;
    end
    
    function Cw = CW(Tw,P)
        % calculates the temperature and concentration at the bubble wall
        thetha = Rv_star/Ra_star*(P./(f_pvsat(Tw*T_inf)/P_inf) -1);
        Cw = 1./(1+thetha);
    end
    
    function Tauw = Boundary(prelim)
        % solves temperature boundary conditions at the bubble wall
        % create finite diff. coeffs.
        % coefficients in terms of forward difference
        
        % second order
        coeff = [-3/2 , 2 ,-1/2 ];
        Tm_trans = Tm(2:3);
        T_trans = flipud(Tau(end-2:end-1));
        
        if masstrans
            % C_trans = flipud(C(end-2:end-1));
            C_trans = flipud(C(end-6:end-1));
            Tauw = chi*(2*iota*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
                chi*(-coeff*[prelim;T_trans])/deltaY + ...
            Fom*L_heat_star*P*( (CW(TW(prelim),P)*(Rv_star-Ra_star)+Ra_star))^-1 *...
                (TW(prelim) * (1-CW(TW(prelim),P))  ).^(-1).*...
            (-coeff*[CW(TW(prelim),P);
            C_trans])/deltaY;
        else
            Tauw = chi*(2*iota*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
                chi*(-coeff*[prelim;T_trans])/deltaY;
        end
    end
    
    function [Diff_Matrix ] = f_finite_diff_mat(Nodes,order,Tm_check)
        % Creates finite difference matrices
        % Nodes: Number of nodes
        % order: order of differentiation ( 1st derivative vs 2nd derivative)
        % Tm_check: 0 not for external temp , 1 used for ext. temp
        
        if Tm_check == 0
            N = Nodes-1;
            deltaY = 1/N;
            K = 1:1:N+1;
            yk = (K-1)*deltaY;
        elseif Tm_check == 1
            N = Nodes-1;
            deltaY = -2/N;
            K = 1:1:N+1;
            yk = 1+(K-1)*deltaY;
        end
        
        Diff_Matrix = zeros(Nodes,1);
        
        if order == 1
            for counter=2:N  %in between
                Diff_Matrix(counter,counter+1) = 1 ;
                Diff_Matrix(counter,counter-1) = -1 ;
            end
            if Tm_check == 0
                Diff_Matrix(end,end) = 3;
                Diff_Matrix(end,end-1) = -4;
                Diff_Matrix(end,end-2) = 1;
            end
            Diff_Matrix = 0.5 * Diff_Matrix / deltaY ;
        elseif order == 2
            if Tm_check == 0
                for counter=2:N  %in between
                    Diff_Matrix(counter,counter+1) = 1 + deltaY/yk(counter);
                    Diff_Matrix(counter,counter)   = -2;
                    Diff_Matrix(counter,counter-1) = 1 - deltaY/yk(counter);
                end
            elseif Tm_check == 1
                for counter=2:N  %in between
                    Diff_Matrix(counter,counter+1) = 1 ;
                    Diff_Matrix(counter,counter)   = -2 ;
                    Diff_Matrix(counter,counter-1) = 1 ;
                end
            end
            if Tm_check == 0
                Diff_Matrix(1,1)= -6;
                Diff_Matrix(1,2) = 6;
            end
            Diff_Matrix = Diff_Matrix / (deltaY^2) ;
        end
        
        Diff_Matrix = sparse(Diff_Matrix);
        
    end
    
end

function cdd = preStressInt(L,N)
    
    cdd = L*N*0;
    
end


% function cdd = StressInt(L,N,varargin)
%
%     if nargin == 2
%         k = 1;
%     else
%         k = varargin{1};
%     end
%     syms x;
%     cdd = zeros(N-k+1,1);
%
%     for n = k:N
%         cdd(n-k+1) = subs(2*L*int((cos(n*acos(x))-1)/((L*(2/(1-x)-1)+1)*(1-x)^2),-1,1));
%     end
%
% end
