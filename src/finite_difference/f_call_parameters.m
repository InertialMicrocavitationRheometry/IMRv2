function [P]  = f_call_parameters(R0,vmaterial)
  % Code to create parameter .mat file for RP_Cav to use 
  % Parameters: 
    A     = 5.28e-5;            % (W/m-K^2)Thermal Conductivity coeff
    B     = 1.165e-2;           % (W/m-K)Thermal Conductivity coeff
    D0    = 24.2e-6;            % Diffusion Coeff m^2/s
    k     = 1.47;               % Ratio of Specific Heats 
    S     = 0.056;              % (N/m) Liquid Surface Tension 
    G     = 1.39E4;             % (Pa) Medium Shear Modulus 
    T_inf = 300;                % (K) Far field temp. 
    [mu8,Dmu,v_a,v_nc,v_lambda] = f_nonNewtonian_Re(vmaterial);
    P_inf = 101325;             % (Pa) Atmospheric Pressure 
    rho   = 999;                %1050.761; % (Kg/m^3) Liquid Density
    Km    = 0.615;              % (W/m-K)Thermal Conductivity Medium
    Cp    = 3.61e3;             % Specific Heat Medium J/Kg K;
    Dm    = Km /(rho*Cp) ;      % Thermal Diffusivity m^2/s 
    Ru    = 8.3144598;          % (J/mol-K) Universal Gas Constant
    Rv    = Ru/(18.01528e-3);   % (J/Kg-K) Gas constant vapor
    Ra    = 438.275;            % (J/Kg-K)Gas constant air
    % Ru/(18.01528e-3);%Ru/(28.966e-3); 
    L     = 1;                  % Strech variable to map domain outside the bubble
    L_heat= 2264.76e3;          % (J/Kg) Latent heat of evaporation 
    C     = 1510;               % sound speed (m/s)
    
  % Intermediate calculated variables
    K_infy  = A*T_inf+B; 
    Rnondim = P_inf/(rho*T_inf);
    Uc      = sqrt(P_inf/rho);
    Pv      = f_pvsat(T_inf );
    P0      = P_inf + 2*S/R0 ;      % need to add Pv_sat at room temp 
    thetha  = Rv/Ra*(P0-Pv)/Pv;     % mass air / mass vapor 
    C0      = 1/(1+thetha); 
    mv0     = Pv*(4/3*pi*R0^3)/Rv/T_inf;
    ma0     = (P0-Pv)*(4/3*pi*R0^3)/Ra/T_inf;
    Mnondim = rho*(4/3*pi*R0^3);
    
  % Final non-dimensional variables
    chi = T_inf*K_infy/(P_inf*R0*Uc);
    fom = D0/(Uc*R0);
    foh = Dm/(Uc*R0); 
    Ca = P_inf/G; 
    Re8 = P_inf*R0/(mu8*Uc);
    if Dmu ~= 0
        DRe = P_inf*R0/(Dmu*Uc);
    else 
        DRe = 0;
    end
    We = P_inf*R0/(2*S);
    Br = Uc^2/(Cp*T_inf);   
    A_star = A*T_inf /  K_infy;
    B_star = B / K_infy; 
    Rv_star = Rv/Rnondim;
    Ra_star = Ra/Rnondim;
    P0_star = P0/P_inf;
    t0 = R0/Uc;
    L_heat_star = L_heat/(Uc)^2;
    Km_star = Km/K_infy; 
    C_star = C/Uc; 
    mv0 = mv0/Mnondim; ma0 = ma0/Mnondim;     
    v_lambda_star = v_lambda/t0;
    
    P = [k chi fom foh Ca Re8 We Br A_star...
         B_star Rv_star Ra_star P0_star t0 C0 L L_heat_star Km_star ...
         P_inf  T_inf C_star mv0  ma0 DRe v_a v_nc v_lambda_star];
end

function [mu8,Dmu,a,nc,lambda] = f_nonNewtonian_Re(vmaterial)
%F_NONNEWTONIAN Outputs the Reynolds number that is dynamically changes
% with the shear rate. Note: Re = P_inf*R0/(m8*Uc). Units are in Pascal
% seconds.
    a = 0; nc = 0; lambda = 0;
    if strcmp('water',vmaterial)==1
        %mu8 = 8.3283e-4;
        %muo = 8.3283e-4;
        mu8 = 0.05;
        muo = 0.05;
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
        error('No viscosity model specified in f_call_parameters, exiting');
    end
    Dmu = muo-mu8;
end