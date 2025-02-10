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
radial          = eqns_opts(1);  bubtherm        = eqns_opts(2); 
medtherm        = eqns_opts(3);  stress          = eqns_opts(4); 
eps3            = eqns_opts(5);  vapor           = eqns_opts(6);
masstrans       = eqns_opts(7);  perturbed       = eqns_opts(8);  
nl              = eqns_opts(9);
if (stress == 4); ptt = 1; else; ptt = 0; end

% solver options
method          = solve_opts(1); spectral        = solve_opts(2); 
divisions       = solve_opts(3); Nv              = solve_opts(4); 
Nt              = solve_opts(5); Mt              = solve_opts(6); 
Lv              = solve_opts(7); Lt              = solve_opts(8); 

% dimensionless initial conditions
Rzero           = init_opts(1);  Uzero           = init_opts(2); 
p0star          = init_opts(3);  P8              = init_opts(4); 
T8              = init_opts(5);  Pv_star         = init_opts(6); 
Req             = init_opts(7);  S0              = init_opts(8);
alphax          = init_opts(9);

% time span options
tspan = tspan_opts;
tfin = tspan(end);
% output options
dimensionalout  = out_opts(1);  progdisplay     = out_opts(2); 

% physical parameters 

% acoustic parameters
Cstar           = acos_opts(1); GAMa            = acos_opts(2); 
kappa           = acos_opts(3); nstate          = acos_opts(4); 
% dimensionless waveform parameters
om              = wave_opts(1); ee              = wave_opts(2); 
tw              = wave_opts(3); dt              = wave_opts(4); 
mn              = wave_opts(5); wave_type       = wave_opts(6); 

pvarargin = [om,ee,tw,dt,mn,wave_type];

% dimensionless viscoelastic
We              = sigma_opts(1); Re8             = sigma_opts(2); 
DRe             = sigma_opts(3); v_a             = sigma_opts(4); 
v_nc            = sigma_opts(5); Ca              = sigma_opts(6); 
LAM             = sigma_opts(7); De              = sigma_opts(8); 
JdotA           = sigma_opts(9); v_lambda_star   = sigma_opts(10); 
iWe             = 1/We;
if Ca==-1; Ca=Inf; end

% dimensionless thermal 
Foh             = thermal_opts(1); Br              = thermal_opts(2); 
alpha           = thermal_opts(3); beta            = thermal_opts(4); 
chi             = thermal_opts(5); iota            = thermal_opts(6); 

% dimensionaless mass transfer 
Fom             = mass_opts(1); C0              = mass_opts(2);
Rv_star         = mass_opts(3); Ra_star         = mass_opts(4);
L_heat_star     = mass_opts(5); mv0             = mass_opts(6);
ma0             = mass_opts(7);

% pre_process

% creates finite difference matrices 
D_Matrix_T_C = f_finite_diff_mat(Nt,1,0);
DD_Matrix_T_C = f_finite_diff_mat(Nt,2,0);
D_Matrix_Tm = f_finite_diff_mat(Mt,1,1);
DD_Matrix_Tm = f_finite_diff_mat(Mt,2,1);

% create spatial nodes

% inside the bubble
N = Nt -1; 
deltaY = 1/N;
i = 1:1:N+1;
y = ((i-1)*deltaY)';
% outside the bubble     
Nm =Mt-1; 
deltaYm = -2/Nm;
j = 1:1:Nm+1;

% precomputations
%LDR = LAM*De/Re8;
sam = 1 - Pv_star + GAMa; 
no = (nstate-1)/nstate;
kapover = (kappa-1)/kappa;
xi = (1+(j-1)*deltaYm)';
yT = ((2./(xi+1)-1)*Lt+1);

% initial condition assembly

% radius, velocity, pressure, bubble temperature, medium temperature,
% vapor concentration
T = zeros(-1,1);
if bubtherm
    Tau0 = zeros(Nt,1);
    Tm0 = ones(Mt ~= -1);
else
    Tau0 = zeros(-1,1);
    Tm0 = zeros(-1,1);
end
if masstrans
    C0 = C0*ones(Nt,1);
else 
    C0 = zeros(-1,1);
end
init = [Rzero; Uzero; p0star; Tau0; Tm0; C0];
tau_del = [];
TL = [];

% solver start
f_display(radial, bubtherm, masstrans, stress, spectral, eps3, Re8, De, Ca, LAM);
stepcount = 0;
bubble = @SVBDODE;
[t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);

% post processing

% extract result
R = X(:,1); 
U = X(:,2); 
p = X(:,3); 
if bubtherm
    Tau = X(:,4:(Nt+3)); 
    T = (A_star -1 + sqrt(1+2*Tau*A_star)) / A_star; % Temp in bubble
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
    stepcount = stepcount + 1;
    if progdisplay == 1, disp(t/tfin); end
    
    % extract standard inputs
    R = X(1); 
    U = X(2); 
    p = X(3); 

    % solve for boundary condition at the wall
    if (medtherm == 1)
        if t/tspan_star> 0.001
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
        Tau = X(4:(Nt+3)); 
        tau_del=[tau_del prelim]; 
        Tau(end) = prelim;
        T = TW(Tau);
         
        K_star = A_star*T+B_star; 
        pVap = (f_pvsat(T(1)*T8)/P8); 
    else
        pVap = p0star;
    end
    if medtherm
        Tm(1) = T(end);
        Tm = X((Mt+4):(2*Mt+3));
    end 

    % updating the viscous forces/Reynolds number    
    % [fnu,intfnu,dintfnu,ddintfnu] = ...
    %     f_nonNewtonian_integrals(vmodel,U,R,v_a,v_nc,v_lambda);

    Taudot = zeros(-1,1);
    Tmdot = zeros(-1,1);
    Cdot = zeros(-1,1);

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
        pdot = 3/R*(chi*(kappa-1)*DTau(end)/R-kappa*p*U+...
              + kappa*p*fom*Rv_star*DC(end)...
              /( T(end)*R*Rmix(end)*(1-C(end)) ) );

        % temperature inside the bubble
        U_vel = (chi/R*(kappa-1)*DTau-yk*R*pdot/3)/(kappa*p);
        first_term = (DDTau*chi/R^2+pdot)*( K_star*T/p*kapover);
        second_term = -DTau*((1/R)*(U_vel-yk*U));
   
        Taudot= first_term+second_term; 
        Taudot(end) = 0;

        % vapor concentration equation 
        U_mix = U_vel + fom/R*((Rv_star - Ra_star)/Rmix)*DC  ;
        one = DDC;
        two = DC*(DTau/(K_star*T)+((Rv_star - Ra_star)/Rmix)*DC );
        three =  (U_mix-U*yk)/R*DC; 

        % concentration evolution
        Cdot = fom/R^2*(one - two) - three;
        Cdot(end) = 0;

    elseif bubtherm 
        % temp. field inside the bubble
        DTau  = D_Matrix_T_C*Tau;
        DDTau = DD_Matrix_T_C*Tau;
        % internal pressure equation
        pdot = 3/R*(chi*(kappa-1)*DTau(end)/R-kappa*p*U);
        % temperature inside the bubble
        U_vel = (chi/R*(kappa-1).*DTau-yk*R*pdot/3)/(kappa*p);
        first_term = (DDTau.*chi./R^2+pdot).*( K_star.*T/p*kapover );
        second_term = -DTau.*(1/(R).*(U_vel-yk*U));
   
        Taudot= first_term+second_term; 
        Taudot(end) = 0;

    else
        % polytropic gas
        p = p0star*(1/R)^(3*kappa);
        pdot= -3*kappa*U/R*p;
        pVap = Pv_star;
    end

    if medtherm
        % temp. field outside the bubble
        DTm = D_Matrix_Tm*Tm;
        DDTm = DD_Matrix_Tm*Tm;
        % warm liquid
        first_term = (1+xi).^2./(Lt*R).*...
            (U./yT.^2.*(1-yT.^3)/2+Foh/R.*((xi+1)/(2*Lt)-1./yT)).* DTm;
        second_term = Foh/R^2.*(xi+1).^4/Lt^2.*DDTm/4;
        %third_term =  4*Br./yT.^6.*(3/Re8.*(U/R)^2);
        third_term =  4*Br./yT.^6.*(3/(Re8+DRe*fnu).*(U/R)^2);
        Tmdot = first_term+second_term+third_term;
        % Sets boundary condition on temp        
        Tmdot(end) = 0; 
        % Previously calculated; 
        Tmdot(1) = 0; 
    end     

    if stress == 0
        % no stress
        J = 0;
        JdotX = 0;
        %TODO Need to add non-Newtonian behavior to JdotX 
        %((1-U/C_star)*R + ...
        %  4/Re8/C_star - 6*ddintfnu*iDRe/C_star);
    elseif stress == 1 
        % Kelvin-Voigt with neo-Hookean elasticity
        J = (4*(Req/R) + (Req/R)^4 - 5)/(2*Ca) - 4/Re8*U/R;
        % JdotX = -2*U*(Req*(1/R)^2 + Req^4*(1/R)^5)/Ca + 4/Re8*U^2/R^2;
    elseif stress == 2
        % quadratic Kelvin-Voigt with neo-Hookean elasticity
        J = (3*alphax-1)*(5 - (Req/R)^4 - 4*(Req/R))/(2*Ca) - 4/Re8*U/R + ...
        (2*alphax/Ca)*(27/40 + (1/8)*(Req/R)^8 + (1/5)*(Req/R)^5 + (1/2)*(Req/R)^2 - ...
        2*R/Req);
        JdotX = ((3*alphax-1)/(2*Ca))*((4*Req^4*U/R^5) + (4*Req*U/R^2)) + ...
        4*(U^2)/(Re8*R^2) - (2*alphax/Ca)*(2*U/Req + Req^8*U/R^9 + ...
        Req^5*U/R^6 + Req^2*U/R^3);       
    else
        % TODO ADD ADDITIONAL STRESS MODELS
        error('stress setting is not available');
    end

    % pressure waveform
    [pf8,pf8dot] = f_pinfinity(t,pvarargin);
    
    % bubble wall acceleration
    [Udot] = f_radial_eq(radial, p, pVap, pf8, pf8dot, iWe, R, U, J, JdotX,...
        Cstar, sam, no, GAMa, nstate, JdotA );

    % stress integral rate
    % TODO ADD STRESS MODELS
    % Jdot = JdotX - JdotA*Udot/R;

    % output assembly
    dXdt = [U; Udot; pdot; Taudot; Tmdot; Cdot];

end
disp('--- COMPLETED SIMULATION ---');    

% functions called by solver 

function Tw= TW(Tauw)
    %calculates the temperature at the bubble wall as a fuction of \tau 
    Tw = (A_star -1 + sqrt(1+2*Tauw*A_star)) / A_star;
end

function Cw= CW(Tw,P)
    % calculates the temperature and concentration at the bubble wall 
    thetha = Rv_star/Ra_star*(P./(f_pvsat(Tw*T_inf)/P_inf) -1);
    Cw = 1./(1+thetha); 
end

function Tauw= Boundary(prelim)       
    % solves temperature boundary conditions at the bubble wall     
    % create finite diff. coeffs. 
    % coefficients in terms of forward difference 
    
    % second order
    coeff = [-3/2 , 2 ,-1/2 ];
    Tm_trans = Tm(2:3);
    T_trans = flipud(Tau(end-2:end-1));
    C_trans = flipud(C(end-2:end-1));
    
    % sixth order    
%     coeff= [-49/20 ,6	,-15/2	,20/3	,-15/4	,6/5	,-1/6]; %Sixth order coeff 
%     Tm_trans = Tm(2:7);
%     T_trans = flipud(Tau(end-6:end-1));
%     C_trans = flipud(C(end-6:end-1));

    Tauw =chi*(2*Km_star/L*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
    chi*(-coeff*[prelim ;T_trans] )/deltaY + Cgrad*...
    fom*L_heat_star*P*( (CW(TW(prelim),P)*(Rv_star-Ra_star)+Ra_star))^-1 *...
    (TW(prelim) * (1-CW(TW(prelim),P))  ).^(-1).*...
    (-coeff*[CW(TW(prelim),P); C_trans] )/deltaY; 
end

end