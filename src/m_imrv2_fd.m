% file m_imrv2_finitediff.m
% brief contains module m_imrv2_finitediff

% brief This module features a fourth- and sixth-order accurate finite
% difference solver of the PDEs involving thermal transport and
% viscoelasticity to solve Rayleigh-Plesset equations
function varargout =  m_imrv2_fd(varargin)
    
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
    Lt               = solve_opts(8);
    
    % dimensionless initial conditions
    Rzero           = init_opts(1);
    Uzero           = init_opts(2);
    pg0star          = init_opts(3);
    P8              = init_opts(4);
    T8              = init_opts(5);
    Pv_star         = init_opts(6);
    Req             = init_opts(7);
    alphax          = init_opts(8);
    R0              = init_opts(9);
    
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
    kapover = (kappa-1)/kappa;
      
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
    ic          = (4+Nt+Mt+Nc):(3+Nt+Mt+Nc+Nv);
    id          = (4+Nt+Mt+Nc+Nv):(3+Nt+Mt+Nc+2*Nv);

    % precomputations for viscous dissipation
    % zT = 1 - 2./(1 + (yT - 1)/Lv);
    cdd = preStressInt(Lv,Nv);

% Need to modify intial conditions for the Out-of-Equilibrium Rayleigh
% Collapse:
% if strcmp(Pext_type,'IC')
    Pv = f_pvsat(1*T8)/P8;
    Pb_star = pg0star + 1*f_pvsat(1*T8)/P8;
    % Need to recalculate intital concentration
    theta = Rv_star/Ra_star*(Pb_star-Pv)/Pv; % mass air / mass vapor
    C0 = 1/(1+theta);
    % Calculate the equilibrium radii ratio for initial stress state:
    % [REq,~,~] = IMRCalc_Req(R0, 1, 1, pg0star*P8, G, G1, mu);

    ma0 = pg0star/Ra_star;
    fun = @(x) Pv*(1+(ma0/x)*(Ra_star/Rv_star))-1-...
    1/We*(Pv/(Rv_star*x))^(1/3) ; % parameterized function

    MTotal0 = pg0star/Ra_star + Pv/Rv_star;
    exp2=MTotal0;
    x = fzero(fun,exp2,optimset('display','off'));

    while (isnan(x))
        exp2 = exp2/1.11;
        x = fzero(fun,exp2,optimset('display','off'));
    end
    MVE = x;
    REq  =(Rv_star*MVE/Pv)^(1/3);
    
    % Zhiren Note: The process of finding REq does not depend on material
    % parameters, [G, G1, mu] are fed mostly to use the IMRcall_Param
    % function ...

    % C0vec = C0*ones(Nt,1);
    %Uzero = -1*(1-Pb_star)/(Cstar); %Intitial velocity due to shockwave
    % Uzero = 0;
    
% end

    % initial condition assembly
    
    % radius, velocity, pressure, bubble temperature, medium temperature,
    % vapor concentration

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
        C0vec = C0*ones(Nc,1);
    else
        C0vec = zeros(-1,1);
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

    tau_del = [];

    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral, ...
        eps3, vapor, Re8, De, Ca, LAM, 'finite difference');
    bubble = @SVBDODE;
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);

    % extract result
    R = X(:,1);
    U = X(:,2);
    p = X(:,3);
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
        U = U*uc;
        p = p*p0;
        if bubtherm
            T = T*T8;
        end
        if medtherm
            Tm = Tm*T8;
        end
    end
    
    % outputs
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = U;
    varargout{4} = p;
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
    function [dXdt] = SVBDODE(t,x)
        
        % extract solution

        % bubble wall radius
        R = x(1); 
        % bubble wall velocity
        U = x(2); 
        % internal bubble pressure
        p = x(3); 
        % auxilary variable for internal bubble temperature
        Tau = x(ibubtherm);
        % temperature in the material
        Tm = x(imedtherm); 
        % vapor concentration inside the bubble
        C = x(imass); 
        % stress integral
        %Z1 = x(4); 

        % if (1 == 1)
            if t/tfin> 0.001
                guess= -.001+tau_del(end);
                prelim  = fzero(@Boundary,guess);
            else
                guess = -.0001;
                prelim  = fzero(@Boundary,guess);
            end
        % else
            % prelim = 0;
        % end
        
        % sets value at boundary conditions
        tau_del = [tau_del prelim];
        Tau(end) = prelim;
        T = TW(Tau);
        Tm(1) = T(end);
        
        % calculate variables
        K_star = alpha*T+beta;
        C(end) =  CW(T(end),p);
        
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

        % internal pressure equation
        pdot = 3/R*(1*chi*(kappa-1)*DTau(end)/R - kappa*p*U +...
              + 1*kappa*p*Fom*Rv_star*DC(end)/( T(end)*R* Rmix(end)* (1-C(end)) ) );
        
        % temperature of the gas inside the bubble
        U_vel = (chi/R*(kappa-1).*DTau-y*R*pdot/3)/(kappa*p);
        first_term = (DDTau.*chi./R^2+pdot).*(K_star.*T/p*kapover);
        second_term = -DTau.*(U_vel-y*U)./R;
        
        Taudot = first_term+second_term;
        Taudot(end) = 0;
        
        % vapor concentration equation
        U_mix = U_vel + Fom/R*((Rv_star - Ra_star)./Rmix).*DC;
        one = DDC;
        two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
        three =  (U_mix-U.*y)/R.*DC;
        
        Cdot = Fom/R^2*(one - two) - three;
        Cdot(end) = 0;
        
        % material temperature equations
        first_term = (1+xi).^2./(Lt*R).*(U./yT.^2.*(1-yT.^3)/2+Foh/R.*((xi+1)/(2*Lt)-1./yT)).* DTm;
        second_term = Foh/R^2.*(xi+1).^4/Lt^2.*DDTm/4;
        third_term =  3*Br./yT.^6.*(4/(3*Ca).*(1-1/R^3)+4.*U/(Re8.*R)).*U./R;
        Tmdot = first_term+second_term+third_term;
        % sets boundary condition on temp
        Tmdot(end) = 0; 
        % previously calculated
        Tmdot(1) = 0; 
        
        % equations of motion
        Rdot = U;
        
        % bubble wall acceleration
        % [Rddot] = f_radial_eq(radial, p, pdot, pVap, pf8, pf8dot, iWe, R, U, ...
            % J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA );

        % pressure waveform
        [pf8,pf8dot] = f_pinfinity(t,pvarargin);

        Pv = (f_pvsat(T(end)*T8)/P8);


        % [J,JdotX,~,~] = ...
            % f_stress_calc(stress,x,Req,R,Ca,De,Re8,U,alphax,ic,id,LAM,zeNO,cdd);
                
                Rst = R/REq;
                J = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re8*U/R ;
                JdotX =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re8*U^2/R^2;
        % if comp == 0
        %     %Rayleigh-Plesset equation
        %     Rddot = (p + abs(1-1)*Pv  - 1 - Pext + J - 1/(We*R) -1.5*U^2)/R;
        % else
            % Keller-Miksis equation
            % if linkv==1 || neoHook==1 || Yeoh==1
                % SdotA = 4/Re8;
            % elseif sls==1 || nhzen==1 || fdkv==1 || zzzen==1 || fdmax==1  % ZZ - Not exactly true for FDKV, but I also don't think this matters ...
                % SdotA = 0;
            % elseif nhkv_pld==1
            %     %SdotA = 4/Re8*(2^alpha+1)/3*(abs(U)/R)^(alpha-1);
            %     SdotA = 4/Re8/3*(2^alpha+1)*sign(U)*(abs(U)/R)^(alpha)*R^2/U^2;
            %     if isnan(SdotA)
            %         SdotA=4/Re8;
            %     end
            % end
            
            % if fdkv == 1
            %     RHS = (1+U/Cstar)...
            %         *(p  + abs(1-1)*Pv -1/(We*R) + J - 1 - Pext)  ...
            %         + R/Cstar*(pdot+ U/(We*R^2) + Jdot -P_ext_prime );
            %     LHS = (3/2)*(1-U/(3*Cstar))*U^2;
            %     denom = (1-U/Cstar)*R - (R/Cstar)*Sdd;
            % 
            %     Rddot = (RHS - LHS)/denom;
            % 
            % else  % Original expression
                 Rddot = ((1+U/Cstar)...
                    *(p  + abs(1-1)*Pv -1/(We*R) + J - 1 - pf8)  ...
                    + R/Cstar*(pdot+ U/(We*R^2) + JdotX - pf8dot) ...
                    - 1.5*(1-U/(3*Cstar))*U^2)/((1-U/Cstar)*R); % +JdotA/(Cstar));
            % end
        % end

        % stress equation
        Z1dot = zeros(-1,1);
        Z2dot = Z1dot;
        % [J,JdotX,Z1dot,Z2dot] = ...
            % f_stress_calc(stress,X,Req,R,Ca,De,Re8,U,alphax,ic,id,LAM,zeNO,cdd);

        Jdot = 0;
        dXdt = [Rdot; 
            Rddot; 
            pdot; 
            Taudot; 
            Tmdot; 
            Cdot;
            Z1dot;
            Z2dot];
         
    end
    % end of solver
    
    % functions called by solver
    
    % calculates the temperature at the bubble wall as a function of \tau
    function Tw = TW(Tauw)
        Tw = (alpha -1 + sqrt(1+2*Tauw*alpha)) / alpha;
    end

    % calculates the concentration at the bubble wall
    function Cw = CW(Tw,p)
        theta = Rv_star/Ra_star*(p./(f_pvsat(Tw*T8)/P8)-1);
        Cw = 1./(1+theta);
    end

    % solves temperature boundary conditions at the bubble wall
    function Tauw = Boundary(prelim)
        
        % coefficients in terms of second-order accurate forward difference
        coeff = [-3/2 , 2 ,-1/2 ];
        Tm_trans = Tm(2:3);
        T_trans = flipud(Tau(end-2:end-1));
        C_trans = flipud(C(end-2:end-1));
                
        Tauw =chi*(2*iota*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
            chi*(-coeff*[prelim ;T_trans] )/deltaY + 1*...
            Fom*L_heat_star*p*( (CW(TW(prelim),p)*(Rv_star-Ra_star)+Ra_star))^-1 *...
            (TW(prelim) * (1-CW(TW(prelim),p))  ).^(-1).*...
            (-coeff*[CW(TW(prelim),p); C_trans] )/deltaY;
        
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
        
        Diff_Matrix = zeros(Nodes);
        
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