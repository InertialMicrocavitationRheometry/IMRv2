function [P]  = IMRcall_parameters(R0,G,G1,mu)

% Code to create parameter .mat file for RP_Cav to use 

% Zhiren - Modified Feb. 2021, add some KVFD parameters
% The idea is that we will spit out multiple parameters. Some will not be
% used anyway.

% Parameters: 

    %0 < G < 10^6; % (Pa) Medium Shear Modulus 
    %0 < mu < 10^3; % (Pa s) Viscocity 
    alpha = 1-eps; % Power law viscosity model exponent 
    G2 = 1E1; %G2 for the Yeoh model
    
    A = 5.28e-5; % (W/m-K^2)Thermal Conductivity coeff
    B = 1.165e-2; % (W/m-K)Thermal Conductivity coeff
    D0 = 24.2e-6; %Diffusion Coeff m^2/s
    k = 1.4; % Ratio of Specific Heats 
    S = 0.072;%56; % (N/m) Liquid Surface Tension 
    T_inf = 298.15; % (K) Far field temp. 
    P_inf = 101325; % (Pa) Atmospheric Pressure 
    rho = 1064;%1060.2; % (Kg/m^3) Liquid Density
    Km = 0.55;  %(W/m-K)Thermal Conductivity Medium
    Cp = 4.181e3; % Specific Heat Medium J/Kg K;
    Dm = Km /(rho*Cp) ; % Thermal Diffusivity m^2/s 
    Ru = 8.3144598; % (J/mol-K) Universal Gas Constant
    Rv = Ru/(18.01528e-3); % (J/Kg-K) Gas constant vapor
    Ra = Ru/(28.966e-3); % (J/Kg-K)Gas constant air
    L = 2; % Strech variable to map domain outside the bubble
    L_heat = 2264.76e3; % (J/Kg) Latent heat of evaporation 
    C = 1484;%1540;%1484; % sound speed (m/s)
    
  % Intermidiate calculated variables
  
    K_infy = A*T_inf+B; 
    Rnondim = P_inf/(rho*T_inf);
    Uc = sqrt(P_inf/rho);
    Pv = Pvsat(T_inf );
    P0 = P_inf + 2*S/R0 ; % need to add Pv_sat at room temp 
    theta = Rv/Ra*(P0-Pv)/Pv; % mass air / mass vapor 
    C0 = 1/(1+theta);

    
  % Final non-dimensional variables

    chi = T_inf*K_infy/(P_inf*R0*Uc);
    fom = D0/(Uc*R0);
    foh = Dm/(Uc*R0); 
    Ca = P_inf/G; 
    Re = P_inf*R0/(mu*Uc);
    We = P_inf*R0/(2*S);
    Br = Uc^2/(Cp*T_inf);
    CaY = P_inf/G2;
    
    A_star = A*T_inf /  K_infy;
    B_star = B / K_infy; 
    Rv_star = Rv/Rnondim;
    Ra_star = Ra/Rnondim;
    P0_star = P0/P_inf;
    t0 = R0/Uc;
    L_heat_star = L_heat/(Uc)^2;
    Km_star = Km/K_infy; 
    C_star = C/Uc; 
    De = (mu/G1)*Uc/R0;
    
    % ZZ - For KVFD: G is still spring, G1 is alpha, while mu = emu = vmu (for
    % now) Intermediate values will also be treated here, just for clarity
    % in future revision ...
    
    FDCa = P_inf/mu;
    FDRe = Re;
    afd = G1;   % Because alpha was used for power-law 
    Ze = (FDCa^(1-afd))*(FDRe^afd); 

    
    P = [k chi fom foh Ca Re We Br A_star...
         B_star Rv_star Ra_star P0_star t0 C0 L L_heat_star Km_star ...
         P_inf  T_inf C_star De alpha CaY ...
         Ze afd]; % ZZ: Add FDKV params.
    
end