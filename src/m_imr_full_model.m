% file m_imrv2_full_model.m
% brief contains module m_imrv2_finitediff

% brief This module features a fourth- and sixth-order accurate finite
% difference solver of the PDEs involving thermal transport and
% viscoelasticity to solve Rayleigh-Plesset equations
function varargout =  m_imr_full_model(varargin)
    
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
    Pb_star         = init_opts(3);
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
    iWe = 1/We;
    if Ca == -1
        Ca = Inf;
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

    % inside the bubble
    N = Nt-1;
    deltaY = 1/N;
    i = 1:1:N+1;
    y = ((i-1)*deltaY)';
    
    % outside the bubble
    Nm = Mt-1;
    deltaYm = -2/Nm;
    j = 1:1:Nm+1;
    xi = (1+(j-1)*deltaYm)';
    yT = ((2./(xi+1)-1)*Lt+1);
    
    % precomputations
    %LDR = LAM*De/Re8;
    sam = 1 - Pv_star + GAMa;
    no = (nstate-1)/nstate;
    chi_kappa = chi / kappa;
    C1_pv = 1.17e11/P8;
    C2_pv = -5200/T8;

    % boundary conditions coefficients
    coeff = [-1.5 , 2 ,-0.5 ];
    Rva_ratio = Rv_star / Ra_star;
    inv_alpha = 1/alpha;
    Rva_diff = Rv_star - Ra_star;
    grad_Tm_coeff = 2*chi*iota/deltaYm*coeff;
    grad_Trans_coeff = -coeff*chi/deltaY;
    grad_C_coeff = -coeff*Fom*L_heat_star/deltaY;    
    
    % index management
    if masstrans == 1
        Nc = Nt;
    else
        Nc = 0;
    end
    if spectral == 0
        Nv = 0;
    end
    if bubtherm == 0
        Nt = 0;
        Mt = 0;
    end
    if bubtherm == 0 && masstrans == 1
        Nt = 0;
    end
    if medtherm == 0
        Mt = 0;
    end
    ibubtherm   = 4:(3+Nt);
    imedtherm   = (4+Nt):(3+Nt+Mt);
    imass       = (4+Nt+Mt):(3+Nt+Mt+Nc);
    ivisco      = (4+Nt+Mt+Nc):(3+Nt+Mt+Nc+Nv);
    ivisco2     = (4+Nt+Mt+Nc+Nv):(3+Nt+Mt+Nc+2*Nv);
    
    % precomputations for viscous dissipation
    % zT = 1 - 2./(1 + (yT - 1)/Lv);
    cdd = preStressInt(Lv,Nv);
    
    % initial condition assembly
        
    % bubble temperature initial condition
    if bubtherm
        Tau0 = zeros(Nt,1);
    else
        Tau0 = [];
    end

    % medium initial temperature
    if medtherm
        Tm0 = ones(Mt,1);
    else
        Tm0 = [];
    end
    
    % mass transfer initial condition
    if masstrans
        C0vec = C0*ones(Nc,1);
    else
        C0vec = [];
    end
    
    % stress spectra
    if stress < 3
        Sp = zeros(2*(Nv - 1)*(spectral == 1),1);
    elseif stress == 3 || stress == 4
        [Sp] = f_max_pre_stress(Req, Cstar, Pv_star, We, Re8, De, Ca, alphax);
    elseif stress == 5
        Sp = zeros(2*(Nv - 1)*(spectral == 1) + 2,1);
    end
    
    % initial condition vector
    init = [Rzero;
    Uzero;
    Pb_star;
    Tau0;
    Tm0;
    C0vec;
    Sp];
    
    guess = -0.0001;

    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral, ...
        eps3, vapor, Re8, De, Ca, LAM, 'finite difference');
    bubble = @SVBDODE;
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);
    
    % extract result
    R    = X(:,1);
    Rdot = X(:,2);
    P    = X(:,3);
    if bubtherm
        Tau = X(:,ibubtherm);
        T = (alpha - 1 + sqrt(1+2*Tau*alpha)) / alpha; % Temp in bubble
        if medtherm
            Tm = X(:,imedtherm);
        end
    end
    
    if masstrans
        C = X(:,imass);
    end
    
    % transform variables back into their dimensional form
    if (dimensionalout == 1)
        t = t*tc;
        R = R*Rref;
        Rdot = Rdot*uc;
        P = P*p0;
        if bubtherm
            T = T*T8;
        end
        if medtherm
            Tm = Tm*T8;
        end
    end
    
    % outputs assembly
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = Rdot;
    varargout{4} = P;
    varargout{5} = T;
    if bubtherm == 1 && medtherm == 1
        varargout{6} = Tm;
    else
        varargout{6} = ((T8 - 1)*dimensionalout + 1)*ones(size(t,1),1);
    end
    if masstrans == 1
        varargout{7} = C;
    end
    
    % solver function
    function [dXdt] = SVBDODE(t,X)
        
        % showing output
        if progdisplay == 1
            disp(t/tfin);
        end
        
        % extract solution
        
        % bubble wall radius
        R       = X(1);
        % bubble wall velocity
        Rdot    = X(2);
        % internal bubble pressure
        P       = X(3);
        % auxiliary variable for internal bubble temperature
        Tau     = X(ibubtherm);
        % temperature in the material
        Tm      = X(imedtherm);
        % vapor concentration inside the bubble
        C       = X(imass);

        prelim = fzero(@Boundary,guess);
        guess = prelim;
        
        % sets value at boundary conditions
        Tau(end) = prelim;
        % using direct function for speed 
        T = inv_alpha*(alpha - 1 + sqrt(1+2*Tau*alpha));
        Tm(1) = T(end);
        
        % calculate variables
        K_star = alpha * T + beta;  
        % using known repeated function here for speed
        C(end) = 1./(1 + Rva_ratio*(P./(C1_pv*exp(C2_pv/T(end))) - 1));
        Rmix = C*Rv_star + (1-C)*Ra_star;
        
        % temperature field of the gas inside the bubble
        DTau  = D_Matrix_T_C*Tau;
        DDTau = DD_Matrix_T_C*Tau;
        
        % concentration of vapor inside the bubble
        DC  = D_Matrix_T_C*C;
        DDC = DD_Matrix_T_C*C;
        
        % temperature field in the material
        DTm = D_Matrix_Tm*Tm;
        DDTm = DD_Matrix_Tm*Tm;
        
        % equations of motion
        
        % internal pressure equation
        Pdot = 3/R*(1*chi*(kappa-1)*DTau(end)/R - kappa*P*Rdot +...
            + 1*kappa*P*Fom*Rv_star*DC(end)/( T(end)*R* Rmix(end)* (1-C(end)) ) );
        
        % temperature of the gas inside the bubble
        U_vel = (chi_kappa/R*(kappa-1).*DTau - y*R*Pdot/3)/(kappa*P);
        first_term = (DDTau.*chi./R^2+Pdot).*(K_star.*T/P*(kappa-1)/kappa);
        second_term = -DTau.*(U_vel-y*Rdot)./R;
        
        Taudot = first_term+second_term;
        Taudot(end) = 0;
        
        % vapor concentration equation
        U_mix = U_vel + Fom/R*((Rv_star - Ra_star)./Rmix).*DC;
        two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
        three =  (U_mix-Rdot.*y)/R.*DC;
        
        Cdot = Fom/R^2*(DDC - two) - three;
        Cdot(end) = 0;
        
        % material temperature equations
        first_term = (1+xi).^2./(Lt*R).*(Rdot./yT.^2.*(1-yT.^3)/2+Foh/R.*((xi+1)/(2*Lt)-1./yT)).* DTm;
        second_term = Foh/R^2.*(xi+1).^4/Lt^2.*DDTm/4;
        third_term =  3*Br./yT.^6.*(4/(3*Ca).*(1-1/R^3)+4.*Rdot/(Re8.*R)).*Rdot./R;
        Tmdot = first_term+second_term+third_term;
        % sets boundary condition on temp
        Tmdot(end) = 0;
        % previously calculated
        Tmdot(1) = 0;
        
        % pressure waveform
        [pf8,pf8dot] = f_pinfinity(t,pvarargin);
        Pv = C1_pv*exp(C2_pv/T(end));
        
        % stress equation
        [J,JdotX,Z1dot,Z2dot] = ...
            f_stress_calc(stress,X,Req,R,Ca,De,Re8,Rdot,alphax,ivisco,ivisco2,LAM,zeNO,cdd);
        
        % bubble wall acceleration
        [Rddot] = f_radial_eq(radial, P, Pdot, (1-vapor)*Pv, pf8, pf8dot, iWe, ...
            R, Rdot, J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA );
        
        % output assembly
        dXdt = [Rdot;
        Rddot;
        Pdot;
        Taudot;
        Tmdot;
        Cdot;
        Z1dot;
        Z2dot];

    end
    % end of solver
    
    % functions called by solver
   
    % solves temperature boundary conditions at the bubble wall
    function Tauw = Boundary(prelim)
        
        Tm_trans = Tm(2:3);
        T_trans = Tau(end-1:-1:end-2);
        C_trans = C(end-1:-1:end-2);
        
        Tw_prelim = inv_alpha*(alpha - 1 + sqrt(1+2*prelim*alpha));
        Pv_prelim = C1_pv*exp(C2_pv/Tw_prelim);
                
        demCw = (1 + (Rva_ratio) * (P ./ Pv_prelim - 1));
        Cw_prelim = 1 / demCw;

        Tauw = grad_Tm_coeff * [Tw_prelim; Tm_trans]  + ...
            grad_Trans_coeff * [prelim ;T_trans] + ...
             grad_C_coeff * P * ((Cw_prelim * Rva_diff + Ra_star))^-1 * ...
            (Tw_prelim * (1 - Cw_prelim))^-1 * [Cw_prelim; C_trans] ;
        
    end
    
    % Creates finite difference matrices
    % Nodes: Number of nodes
    % order: order of differentiation ( 1st derivative vs 2nd derivative)
    % Tm_check: 0 not for external temp , 1 used for external temp
    function [Diff_Matrix] = f_finite_diff_mat(Nodes,order,Tm_check)
        
        % coordinate creation
        if Tm_check == 0
            N = Nodes-1;
            deltaY = 1/N;
            K = 1:1:N+1;
            y = (K-1)*deltaY;
        elseif Tm_check == 1
            N = Nodes-1;
            deltaY = -2/N;
            K = 1:1:N+1;
            y = 1+(K-1)*deltaY;
        end
        
        Diff_Matrix = zeros(Nodes, Nodes);
        
        if order == 1
            % in between
            for counter = 2:N
                Diff_Matrix(counter,counter+1) = 0.5 ;
                Diff_Matrix(counter,counter-1) = -0.5 ;
            end
            if Tm_check == 0
                Diff_Matrix(end,end) = 1.5;
                Diff_Matrix(end,end-1) = -2;
                Diff_Matrix(end,end-2) = 0.5;
            end
            Diff_Matrix = Diff_Matrix / deltaY ;
            
        elseif order == 2
            if Tm_check == 0
                % in between
                for counter = 2:N
                    Diff_Matrix(counter,counter+1) = 1 + deltaY/y(counter);
                    Diff_Matrix(counter,counter)   = -2;
                    Diff_Matrix(counter,counter-1) = 1 - deltaY/y(counter);
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
    
    function cdd = preStressInt(Lt,N)
        cdd = Lt*N*0;
    end
    
end
