function [ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] =...
    m_cavitation(tspan,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,...
            Tgrad,Tmgrad,Cgrad,Dim,comp,vmaterial,vmodel)

% Original Authors: Carlos Barajas
% Developer(s): Mauro Rodriguez (mauro_rodriguez@brown.edu)
%
% Description: This code is a reduced version of the IMR code taken from
% Estrada et al. (2018) JMPS. Additional physics have been added including
% the Keller-Miksis with enthalpy and non-Newtonian viscosity. 

% Inputs: 
% tspan - time to run simulation
% R0 - Initial Radii
% NT - number of nodes for temperature and concentration fields 
% NTM - number of nodes for temperature in the medium 
% Pext_type - type of external pressure ('sn' = sine, 'RC' = Rayleigh 
% collapse, 'RG' = Rayleigh growth, impulse 'ip',...
% non-equlibrium initial conditions (laser caviation and Flynn(1975) ) 'IC'
% Pext_Amp_Freq - amplitude and frequency of external pressure [amp w]

% Note: For this code the out-of-equilibrium Rayleigh Collapse the intial
% mass in the bubble and radii are specified 


% FOR THE FOLLOWING INPUTS 0 = FALSE AND 1 = TRUE 
% disptime - Displays elapsed time on the command window
% Tgrad - Models temperature gradients in the buuble
% Tmgrad- Models temperature gradients outside the buuble
% Cgrad - Models concentration gradients in the buuble

% Outputs: 
% t - time vector
% T_Bubble - Temperature inside the bubble  
% T_Medium - Temperature outside the bubble  
% R - Bubble Radius 
% U - Bubble velocity 
% P - Internal bubble pressure
% C - Vapor Concentration in the bubble
% Tm - Temperature in the medium 
% Dim - outputs variables in dimensional form
% Comp - 0 (ignores compressibility effects) or 1 (uses Keller- Miksis)
% Reduced - 0 utilizes full model or 1 uses Preston's reduced order model

%***************************************
% Load Parameters : 
Pmt = f_call_parameters(R0,vmaterial); % Calls parameters script 
k = Pmt(1); chi = Pmt(2);  fom = Pmt(3); foh = Pmt(4);  Ca = Pmt(5);  
Re8 = Pmt(6); We = Pmt(7);  Br = Pmt(8);  A_star = Pmt(9); B_star = Pmt(10);
Rv_star = Pmt(11);  Ra_star = Pmt(12); P0_star = Pmt(13); t0 = Pmt(14);
C0 = Pmt(15); L = Pmt(16); L_heat_star = Pmt(17); Km_star = Pmt(18); 
P_inf = Pmt(19); T_inf = Pmt(20); C_star = Pmt(21);  
Mv0 = Pmt(22);   Ma0 = Pmt(23); 
% additional non-Newtonian parameters
DRe = Pmt(24); v_a = Pmt(25); v_nc = Pmt(26); v_lambda = Pmt(27);
if DRe==0
    iDRe = 0;
else
    iDRe = 1/DRe;
end

%****************************************
% Needed to account for fast diffusion

    P0_star = P0_star - (1-Cgrad)*f_pvsat(1*T_inf)/P_inf; 
    % When we assume water vapor undergoes infinitely fast mass diffusion
    % the vapor pressure is constant and P is the pressure of
    % non-condesible gas 

%******************************************
% Creates finite difference matrices 
D_Matrix_T_C = f_finite_diff_mat(NT,1,0);
DD_Matrix_T_C = f_finite_diff_mat(NT,2,0);
D_Matrix_Tm = f_finite_diff_mat(NTM,1,1);
DD_Matrix_Tm = f_finite_diff_mat(NTM,2,1);

%******************************************
% Create spatial nodes
% Inside the bubble
N = NT -1; 
deltaY = 1/N;
i = 1:1:N+1;
yk = ((i-1)*deltaY)';
% Outside the bubble     
Nm =NTM-1; 
deltaYm = -2/Nm;
j = 1:1:Nm+1;
xk = (1+(j-1)*deltaYm)';
yk2 = ((2./(xk+1)-1)*L+1);

%******************************************
% Initial Conditions
tspan_star = tspan/t0;
tspan_star = 5; % forcing this for 3dasm project
R0_star = 1; 
U0_star = 0;  % Change as needed 
Tau0 = zeros(1,NT);
C0 = C0*ones(1,NT);
Tm0 = ones(1,NTM);

% Need to modify initial conditions for the Out-of-Equilibrium Rayleigh
% Collapse:  
if  (Pext_type == 'IC') 
    Pv = f_pvsat(1*T_inf)/P_inf;
    P0_star = Pext_Amp_Freq(1)/P_inf + Cgrad*f_pvsat(1*T_inf)/P_inf; 
    % Need to recalculate intital concentration
    thetha = Rv_star/Ra_star*(P0_star-Pv)/Pv; % masp air / mass vapor 
    C0 = 1/(1+thetha);
    Ma0 = (P0_star-Pv)/Ra_star;  
    % Need to calculate the equilibrium radii for initial 
    % stress state: 
    [REq] = f_calc_Req(R0, Tgrad ,Cgrad,Pext_Amp_Freq(1),vmaterial);
    %REq = 1;
    C0 = C0*ones(1,NT);
    U0_star = -(1-P0_star)/(C_star); % Initial velocity 
    %Plesset & Prosperetti, ARFM 1977, p166
    else
    U0_star = 0;
end

if  (Pext_type == 'RC') 
    U0_star = -1*(Pext_Amp_Freq(1)/P_inf)/(C_star); % Intitial velocity 
     %Plesset & Prosperetti, ARFM 1977, p166 
end

tau_del= [];
tdel=[];
Tdel = []; 
Cdel = []; 

if (Pext_type == 'HN')
    in = load('./data/workspace.mat','pp_HN');
    pp_HN = in.pp_HN;
    in = load('./data/workspace.mat','dp_HN');
    dp_HN = in.dp_HN;
elseif (Pext_type == 'MN')
    in = load('./data/workspace.mat','pp_MN');
    pp_MN = in.pp_MN;
    in = load('./data/workspace.mat','dp_MN');
    dp_MN = in.dp_MN;
elseif (Pext_type == 'ML')
    in = load('./data/workspace.mat','pp_ML');
    pp_ML = in.pp_ML;
    in = load('./data/workspace.mat','dp_ML');
    dp_ML = in.dp_ML;
else
end

%************************************************
% March equations in time 
% options = odeset('RelTol',1e-10); %Tune tolerances of simulation
% [t,X] = ode45(@fun,time_span,X0,options) <== Syntax 
    opts = odeset('RelTol',1e-8,'AbsTol',1E-8);
    X0 = [R0_star U0_star P0_star Tau0 C0 Tm0 ]; 
    [t , X] = ode23tb(@bubble, [0 tspan_star] , X0, opts);
    %[t , X] = ode23tb(@bubble, [0 tspan_star] , X0);
    R = X(:,1); % Bubble wall Radius 
    U = X(:,2); % Bubble wall velocity
    P = X(:,3); % Internal pressure
    Tau = X(:,4:(NT+3)); % Variable relating to internal temp
    C =  X(:,(NT+4):(2*NT+3)); % Vapor concentration in the bubble 
    Tm = X(:, (2*NT+4):end ); % Temperature variation in the medium
    T = (A_star -1 + sqrt(1+2*Tau*A_star)) / A_star; % Temp in bubble
    Mv = 0;

% ******************************
% Transform variables back into their dimensional form 
 if (Dim == 1)
    R = R*R0; 
    t= t*t0;
    T = T*T_inf;
    P = P*P_inf;
    U = U*(R0/t0);
    tdel= tdel*t0; %time vector
    Tdel = Tdel*T_inf;
    Mv = Mv*(Mv0+Ma0);
 end
 %***********************

%*************************************************************************
% Nested function; ODE Solver calls to march governing equations in time
% This function has acess to all parameters above; 

% Solves the full model (PDE's)
function dxdt = bubble(t,x)
    
     % Break x vector into indv. values
     R = x(1); % Bubble wall Radius 
     U = x(2); % Bubble wall velocity
     P = x(3); % Internal pressure
     Tau = x(4:(NT+3));
     C = x((NT+4):(2*NT+3)); 
     Tm = x((2*NT+4):end );
     
    if (disptime == 1 )
        %display simulation time in terminal 
        disp(t/tspan_star);
    end
         
    % *********Solves for boundary condition at the wall************** 
    
    if (Tmgrad == 1)
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

    %****************************************************************
     % Sets value at boundary conditions
     tau_del=[tau_del prelim]; 
     Tau(end) = prelim;
     T = TW(Tau);
     Tm(1) = T(end); 
           
     % Calculated variables     
     K_star = A_star*T+B_star;  
     C(end) =  CW(T(end),P);
     
     Rmix = C*Rv_star + (1-C)*Ra_star; 
     
     % Gets variables that are not directly calculated as outputs
     Tdel = [Tdel , T(end)];
     tdel = [tdel t];
     Cdel = [Cdel C(end)];    
     
     %Set external pressure
     if (Pext_type == 'sn')
         Pext =  -Pext_Amp_Freq(1)/P_inf*sin(2*pi*Pext_Amp_Freq(2)*t*t0) ; 
         P_ext_prime = -2*pi*Pext_Amp_Freq(2)*t0*Pext_Amp_Freq(1)/P_inf...
             *cos(2*pi*Pext_Amp_Freq(2)*t*t0) ;
      
     elseif (Pext_type == 'RC')
         
         Pext = Pext_Amp_Freq(1)/P_inf ; 
                P_ext_prime = 0; 
                
     elseif (Pext_type == 'GS')
         a = 10;
         tw = 1;
         Pext = -Pext_Amp_Freq(1)*(exp(-((t-a)/tw)^2))/P_inf;
         P_ext_prime = 2*Pext_Amp_Freq(1)/P_inf.*(exp(-((t-a)./tw).^2)).*(t*t0-a)./tw^2;

     elseif (Pext_type == 'RG')
         
         Pext = -Pext_Amp_Freq(1)/P_inf ;              
         P_ext_prime = 0;  
         
     elseif (Pext_type == 'ip')
         
         Pext = -Pext_Amp_Freq(1)/P_inf*... 
         (1-heaviside(t-Pext_Amp_Freq(2)/t0)) ;
         P_ext_prime = 0; 
         
      elseif (Pext_type == 'IC')
         
         Pext = 0;
         P_ext_prime = 0; 
         
     elseif (Pext_type == 'HN')
         Pext = ppval(pp_HN,t);
         P_ext_prime = ppval(dp_HN,t);
         
     elseif (Pext_type == 'MN')
         Pext = ppval(pp_MN,t);
         P_ext_prime = ppval(dp_MN,t);
         
     elseif (Pext_type == 'ML')
         Pext = ppval(pp_ML,t);
         P_ext_prime = ppval(dp_ML,t);
         
     end

    % *****************************************    
    % Create derivative terms
    
    % Temp. field inside the bubble
    DTau  = D_Matrix_T_C*Tau;
    DDTau = DD_Matrix_T_C*Tau;
    
    % Concentration field inside the bubble
    DC  = D_Matrix_T_C*C; 
    DDC = DD_Matrix_T_C*C;
    
    % Temp. field outside the bubble
 
    DTm = D_Matrix_Tm*Tm;
    DDTm = DD_Matrix_Tm*Tm;
 
    %***************************************
    % Internal pressure equation
    
     P_prime = 3/R*(Tgrad*chi*(k-1)*DTau(end)/R-k*P*U+...
              + Cgrad*k*P*fom*Rv_star*DC(end)...
              /( T(end)*R* Rmix(end)* (1-C(end)) ) );

    % *****************************************

    %***************************************
    % Updating the viscous forces/Reynolds number    
     [fnu,intfnu,dintfnu,ddintfnu] = ...
         f_nonNewtonian_integrals(vmodel,U,R,v_a,v_nc,v_lambda);
    
    %***************************************    
    % Temperature inside the bubble
    
    U_vel = (chi/R*(k-1).*DTau-yk*R*P_prime/3)/(k*P);
    first_term = (DDTau.*chi./R^2+P_prime).*( K_star.*T/P*(k-1)/k);
    second_term = -DTau.*(1/(R).*(U_vel-yk*U));
   
    Tau_prime = first_term+second_term; 
    Tau_prime(end) = 0;
    Tau_prime = Tau_prime*Tgrad; 
    % *****************************************
      
    %***************************************
    % Vapor concentration equation 
    
    U_mix = U_vel + fom/R*((Rv_star - Ra_star)./Rmix).*DC  ;
    one = DDC;
    two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
    three =  (U_mix-U.*yk)/R.*DC; 
    % *****************************************
    
    % *****************************************
    % C_prime
    
    C_prime = fom/R^2*(one - two) - three;
    C_prime(end) = 0;
    C_prime = C_prime*Cgrad;
    %***************************************
    
    %***************************************
    % External temperature: In the liquid
    
    first_term = (1+xk).^2./(L*R).*...
        (U./yk2.^2.*(1-yk2.^3)/2+foh/R.*((xk+1)/(2*L)-1./yk2)).* DTm;
    second_term = foh/R^2.*(xk+1).^4/L^2.*DDTm/4;
    %third_term =  4*Br./yk2.^6.*(3/Re8.*(U/R)^2);
    third_term =  4*Br./yk2.^6.*(3/(Re8+DRe*fnu).*(U/R)^2);
    Tm_prime = first_term+second_term+third_term;
    Tm_prime(end) = 0; % Sets boundary condition on temp        
    Tm_prime(1) = 0; % Previously calculated; 
    Tm_prime = Tm_prime*Tgrad;
    % ***************************************** 

    %***************************************
    % Viscoelastic Forces : 
    
    E = 1/(2*Ca)*(5-4/R - 1/R^4)+ 4*U/(Re8*R) - 6*intfnu*iDRe; 
    E_prime = 2*U*(1/R^2 + 1/R^5)/Ca - 4/Re8*U^2/R^2 - 6*dintfnu*iDRe; 
     
    %     % Yang and Church (2005); linear elasticity 
    %      E = 4/(3*Ca)*(1 - 1/R^3) + 4/Re*U/R;
    %      E_prime= 4/Ca*U*/R^4 - 4/Re*U^2/R^2;
     
     if  (Pext_type == 'IC') 
     %Account for initial stress state     
       E = 1/(2*Ca)*(5-(4*REq/R) - (REq/R)^4)+ 4*U/(Re8*R); 
       E_prime = 2*U*(REq/R^2 + REq^4/R^5)/Ca - 4/Re8*U^2/R^2;
     
    %     % Yang and Church (2005); linear elasticity 
    %      E = 4/(3*Ca)*(1 - REq^3/R^3) + 4/Re*U/R;
    %      E_prime= 4/Ca*U*REq^3/R^4 - 4/Re*U^2/R^2;
    %      
     end
    %****************************************************

    % Equations of motion 
     R_prime = U;  
  
     if (Tgrad == 0)
        P = P0_star*(1/R)^(3*k);
        P_prime = -3*k*U/R*P;
     end

     Pv = (f_pvsat(T(end)*T_inf)/P_inf);
     if comp == 0
     %R.P
        U_prime = (P +abs(1-Cgrad)*Pv-1-Pext - E - 1/(We*R) - 1.5*U^2)/R;
	 else
     % Keller-Miksis in pressure
        U_prime = ((1+U/C_star)...
          *(P  + abs(1-Cgrad)*Pv - 1 - Pext - E - 1/(We*R)) ...
          + R/C_star*(P_prime+ U/(We*R^2) - E_prime - P_ext_prime ) ...
          - 1.5*(1-U/(3*C_star))*U^2)/((1-U/C_star)*R + ...
          4/Re8/C_star - 6*ddintfnu*iDRe/C_star);
     end
     % *****************************************
     
     dxdt = [R_prime; U_prime; P_prime; Tau_prime; C_prime; Tm_prime]; 
end

%*************************************************************************

% Other nested functions used to carry out repeatetive calculations 
% throughout the code 
function Tw= TW(Tauw)
 
    %calculates the temperature at the bubble wall as a fuction of \tau 
    
    Tw = (A_star -1 + sqrt(1+2*Tauw*A_star)) / A_star;
    
end

function Cw= CW(Tw,P)
    
  % Calculates the concentration at the bubble wall 
  %Function of P and temp at the wall 
 
  thetha = Rv_star/Ra_star*(P./(f_pvsat(Tw*T_inf)/P_inf) -1);
  Cw = 1./(1+thetha); 
  
end

function Tauw= Boundary(prelim)
       
    % Solves temperature boundary conditions at the bubble wall     
    % Create finite diff. coeffs. 
    % Coefficients in terms of forward difference 
    
%    %Second order
    coeff = [-3/2 , 2 ,-1/2 ];
    Tm_trans = Tm(2:3);
    T_trans = flipud(Tau(end-2:end-1));
    C_trans = flipud(C(end-2:end-1));
    
%   Could implement any order... sixth order example is shown below 

%     Sixth order    
%     coeff= [-49/20 ,6	,-15/2	,20/3	,-15/4	,6/5	,-1/6]; %Sixth order coeff 
%     Tm_trans = Tm(2:7);
%     T_trans = flipud(Tau(end-6:end-1));
%     C_trans = flipud(C(end-6:end-1));

    Tauw =chi*(2*Km_star/L*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
    chi*(-coeff*[prelim ;T_trans] )/deltaY + Cgrad*...
    fom*L_heat_star*P*( (CW(TW(prelim),P)*(Rv_star-Ra_star)+Ra_star))^-1 *...
    (TW(prelim) * (1-CW(TW(prelim),P))  ).^(-1).*...
    (-coeff*[CW(TW(prelim),P); C_trans] )/deltaY; 
 %************************************************************************
end
end
