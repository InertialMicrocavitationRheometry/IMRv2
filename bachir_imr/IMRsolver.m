function [t,R,U,P,S,T,C,Tm,tdel,Tdel,Cdel] = IMRsolver(model,G,G1,mu,tspan,R0,NT,NTM, ...
    Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,IX,RMesh)

% Authors:
% Carlos Barajas
% carlobar@umich.edu
% Umich Mechanical Engineering BS '16
% Jon Estrada
% jonathan_estrada@brown.edu
% Brown Solid Mechanics, PhD '17

% Last JBE Update: 8/21/2017


% EDIT - 2021
%
% Zhiren Zhu, zhiren@umich.edu
%
% Added Features:
% 1) Inclusion of FDKV model: Added input IX for initial expansion, FD
% solvers, etc.
% 2) Use optimized spatial mesh for FDKV model - Feb. 2021
% 3) Re-arrange solver code structure, transfer nested functions to local
% functions (pass variables properly, without sharing global variables ...)
% 4) Flipped Astore to allow let spatial node be row, time step be column +
% pre-allocate to prevent change of dimension during solution process.
% 5) Modify Zener model to consider stress in Maxwell element dye to
% initial expansion. The IX passed over for FD model is utilized. - Apr. 2021
% 6) Revised Zener model added as an extra option. Use numerical spatial
% integration to evaluate stress integral. - May. 2021
% 7) Add Maxwell approximation of FD (mainly in solver) - Oct. 2021 (Still
% under development, as of Dec. 2021)
% 8) (Ongoing - Jan.'22) Downsize BT structure for fdmax, so that storage
% is not so heavy? Possible?

% Inputs:
% tspan - time to run simulation
% R0 - Initial Radii
% NT - number of nodes for temperature and concentration fields
% NTM - number of nodes for temperature in the medium
% Pext_type - type of external pressure ('sn' = sine, 'RC' = Rayleigh
% collapse, 'RG' = Rayleigh growth, impulse 'ip',...
% non-equlibrium initial conditions (laser caviation and Flynn(1975) ) 'IC'
% Pext_Amp_Freq - amplitude and frequency of external pressure [amp w]
% IX - [Nx4]: contains time, radius, velocity, acceleration for initial
% expansion, to be used for calculation of FD in FDKV model.

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

%********************************************************************
% Citations for code:
%M.T. Warnez and E. Johnsen, "Numerical modeling of bubble dynamics in
%viscoelastic media with relaxation," Phys. Fluids 27, (2015).

%R. Gaudron, M.T. Warnez, and E. Johnsen, "Bubble dynamics in a
%viscoelastic medium with nonlinear elasticity,"
%J. Fluid Mech. 766, 54-75 (2015).

%X. Yang and C.C. Church, "A model for the dynamics of gas bubbles
%in soft tissue," J.Acoust. Soc. Am. 118, 3595-3606 (2005).

%A. Prosperetti, L. A. Crum, and K.W. Commander, "Nonlinear bubble
%dynamics," J.Acoust. Soc. Am. 83, 502-514 (1988).

%A.T. Preston, "Modeling heat and mass transfer in bubbly cavitating
%flows and shockwaves in cavitating nozzles," Ph.D. thesis,
%California Institute of Technology (2004).

%R.I. Nigmatulin, N.S. Khabeev, and F.B. Nagiev, "Dynamics, heat and mass
%transfer of vapour-gas plobubbles in a liquid," Int. J.
%Heat Mass Transfer, 24, 1033-1044 (1981).
%*************************************************************************

%***************************************
% Load Parameters :
Pmt = IMRcall_parameters(R0,G,G1,mu); % Calls parameters script
k = Pmt(1); 
chi = Pmt(2);
fom = Pmt(3);
foh = Pmt(4); 
Ca = Pmt(5);
Re = Pmt(6);
We = Pmt(7);
Br = Pmt(8); 
A_star = Pmt(9);
B_star = Pmt(10);
Rv_star = Pmt(11); 
Ra_star = Pmt(12);
P0_star = Pmt(13);
t0 = Pmt(14);
C0 = Pmt(15);
L = Pmt(16);
L_heat_star = Pmt(17);
Km_star = Pmt(18);
P_inf = Pmt(19);
T_inf = Pmt(20);
C_star = Pmt(21);
De = Pmt(22);
alpha = Pmt(23);
CaY = Pmt(24);
Ze = Pmt(25);
afd = Pmt(26);


% Zhiren Notes: Ze is "Zener #" = (Ca^(1-afd))*(Re^(1-afd)) in FD springpot
%****************************************

% Material Choice
neoHook = 0;
nhzen = 0;
sls = 0;
linkv = 0;
nhkv_pld = 0;
Yeoh = 0;
fdkv = 0;
zzzen = 0;
fdmax = 0;

if strcmp(model,'neoHook') == 1
    neoHook = 1;
elseif strcmp(model,'Yeoh') == 1
    Yeoh = 1;
elseif strcmp(model,'nhzen') == 1
    nhzen = 1;
elseif strcmp(model,'sls') == 1
    sls = 1;
elseif strcmp(model,'linkv') == 1
    linkv = 1;
elseif strcmp(model,'nhkv_pld') == 1
    nhkv_pld = 1;
elseif strcmp(model,'zzzen') == 1
    zzzen = 1;    
elseif strcmp(model,'fdkv') == 1
    fdkv = 1;
else
    fdmax = 1;
end


% Needed to account for fast diffusion
P0_star = P0_star - (1-Cgrad)*Pvsat(1*T_inf)/P_inf;

% When we assume water vapor undergoes infinitely fast mass diffusion
% the vapor pressure is constant and P is the pressure of
% non-condesible gas

%******************************************
% Creates finite difference matrices
D_Matrix_T_C = Finite_diff_mat(NT,1,0);
DD_Matrix_T_C = Finite_diff_mat(NT,2,0);
D_Matrix_Tm = Finite_diff_mat(NTM,1,1);
DD_Matrix_Tm = Finite_diff_mat(NTM,2,1);
%******************************************

%******************************************
% Create spatial nodes

% Inside the bubble
N = NT-1;
deltaY = 1/N;
i = 1:1:N+1;
yk = ((i-1)*deltaY)';

% Outside the bubble
Nm = NTM-1;
deltaYm = -2/Nm;
j = 1:1:Nm+1;
xk = (1+(j-1)*deltaYm)';
yk2 = ((2./(xk+1)-1)*L+1);

%******************************************
% Initial Conditions
tspan_star = tspan/t0;
R0_star = 1;
U0_star = 0;  % Change as needed
%Z10 = 0;
S0 = 0;
Tau0 = zeros(1,NT);
C0 = C0*ones(1,NT);
Tm0 = ones(1,NTM);
if strcmp(Pext_type,'ga')
    dt_star = Pext_Amp_Freq(2)/t0;
    tw_star = Pext_Amp_Freq(3)/t0;
end

% Need to modify intial conditions for the Out-of-Equilibrium Rayleigh
% Collapse:
if strcmp(Pext_type,'IC')
    Pv = Pvsat(1*T_inf)/P_inf;
    P0_star = Pext_Amp_Freq(1)/101325 + Cgrad*Pvsat(1*T_inf)/P_inf;
    % Need to recalculate intital concentration
    theta = Rv_star/Ra_star*(P0_star-Pv)/Pv; % mass air / mass vapor
    C0 = 1/(1+theta);
    
    % Calculate the equilibrium radii ratio for initial stress state:
    [REq,~,~] = IMRCalc_Req(R0, Tgrad, Cgrad, Pext_Amp_Freq(1), G, G1, mu);
    %REq = 1; %removed 6/15/16 by Jon
    
    % Zhiren Note: The process of finding REq does not depend on material
    % parameters, [G, G1, mu] are fed mostly to use the IMRcall_Param
    % function ...
    
    C0 = C0*ones(1,NT);
    %U0_star = -1*(1-P0_star)/(C_star); %Intitial velocity due to shockwave
    U0_star = 0;
    if sls == 1 || linkv == 1
        S0 = -4/(3*Ca)*(1-REq^3);
    elseif neoHook == 1 || nhkv_pld == 1 % || nhzen == 1  -- Taken out by Zhiren
        S0 = -1/(2*Ca)*(5-REq^4-4*REq);
        
    elseif nhzen == 1
        % Run maxwell_expand to get elastic stretch at end of expansion
        
        % Operate with absolute scale first:
        lam0 = 1/REq; % Stretch at r = R, t = 0;
        texp = -IX(1,1); % Effective duration of expansion
        tau = mu/G1; % Inherent time scale in Maxwell element;
        
        % We need to perform spatial integration to get stress integral in
        % Maxwell element. 
        
        % Borrow spatial mesh from FD model
        RR = 1 + RMesh; % Note: mesh stored as (r0-1) in ref. config.
        RR = RR*REq; % Multiply everything by actual bubble radius
        nr = length(RR);
        
        [SR,RT] = maxwell_stress(RR,lam0,texp,tau); % Ignore modulus for now
        
        Sneq = (1/De)*(1/Re)*trapz(RT,SR);
        
        % Then add ground-state elastic part to get total:
        S0 = -1/(2*Ca)*(5-REq^4-4*REq) + Sneq;
        
    elseif fdkv == 1
        
        % First, we need to scale all the IX params:
        Uc = R0/t0; % Material speed. I think this is convenient to have.
        Ac = Uc/t0; % Characteristic acceleration. Again, just for convenience.
        
        TX = IX(:,1)/t0;
        RX = IX(:,2)/R0;
        VX = IX(:,3)/Uc;
        AX = IX(:,4)/Ac;
        
        % Next, generate spatial mesh, since we need integrate spatially ...
        % Use 3-domain setup to ensure good accuracy:
        
        % Feb. 2021 - Retire old mesh
%         nr1 = 1500;
%         rmax1 = 4.3E-5; 
%         mesh1 = 0:rmax1/(nr1-1):rmax1; 
%         RR1 = 10.^mesh1;
%         nr2 = 1500; 
%         rmax2 = 0.042;
%         mesh2 = rmax1:(rmax2-rmax1)/nr2:rmax2;
%         RR2 = 10.^mesh2(2:end); % Drop first term, repetitive
%         nr3 = 200; 
%         rmax3 = 2;
%         mesh3 = rmax2:(rmax3-rmax2)/nr3:rmax3;
%         RR3 = 10.^mesh3(2:end); % Drop first term, repetitive 
%         RR = [RR1,RR2,RR3]; % This is actually R/R0 (1 to max)
        
        % Passed from outside:
        RR = 1 + RMesh; % Note: mesh stored as (r0-1) in ref. config.

        RR = RR*REq; % Multiply everything by actual bubble radius
        nr = length(RR);
        
        % Set up history storage of FD related terms:
        nchunk = 5000; % Some big size to pre-allocate, updating size of 3D matrix is not fun.
        GT = zeros(nchunk,nr,8);
        
        npre = length(TX); % # of steps in initial expansion history
        
        Astore = zeros(nr,nchunk); % zeros(npre,nr); % For storage of acceleration (Flipped + expand 2/25/21)
        
        for kk = 1:npre
            rbk = RX(kk);
            vbk = VX(kk);
            abk = AX(kk);
            
            % Get all deformed radius in matrix operation:
            RK = ((rbk^3-REq^3) + RR.^3).^(1/3);
            
            % Get current velocity accordingly:
            VK = vbk*(rbk^2)*(RK).^(-2);
            
            % Get acceleration:
            RRBJ = rbk./RK; 
            Astore(:,kk) = abk*(RRBJ.^2) + 2*(VK.^2).*RRBJ./RK.*(1-RRBJ.^3); 
            
            % Go through spatial nodes
            for jj = 1:nr
                % (2/25/21 - converted to matrix operations)
                %rj0 = RR(jj);   % Undeformed radius at node
                %rjk = (rbk^3+rj0^3-REq^3)^(1/3); % Current radius at node
                %vjk = vbk*(rbk/rjk)^2;  % Current velocity at node
                
                rjk = RK(jj);
                vjk = VK(jj);
                
                % Get FD related values: (Using embedded function below)
                GJ = evalgt(rjk,vjk);
                
                % Update storage:
                for ii = 1:8
                    GT(kk,jj,ii) = GJ(ii);
                end % i - FD-related terms
                
                % Additionally, update acceleration: 
                % (2/25/21 - converted to matrix operations; flipped Astore)
                % rrbj = rbk/rjk; % Current radii ratio between bubble wall & node
                % Astore(kk,jj) = abk*rrbj^2+2*(vbk^2)*rrbj/rjk*(1-rrbj^3); 
            end % j - spatial nodes
        end % k - temporal steps
        
                
        Tstep = TX; % Initialize for all FD cases
        
        % Now, get initial stress
        if afd == 0
            S0 = -(1/(2*Ca)+1/(2*Ze))*(5-REq^4-4*REq); % Ze is the Ca of springpot in this case
            GTnow = []; % Generate it to pass on
        else
            % Initialize GTnow, storing FD-related values for current configuration
            GTnow = zeros(8,nr);
            for jj = 1:nr
                for ii = 1:8
                    GTnow(ii,jj) = GT(npre,jj,ii);
                end
            end
            % GTnow will be updated as we march ODE forward
            
            % Define R & U here:
            R = 1;
            U = 0;
            
            % Evaluate stress integral: recommend using "caputo_exact".
            [SJ,~,~,KR] = caputo_exact(Tstep,GT,GTnow,0,RR, REq, R, U, Astore, afd); % set t=0
            %[SJ,dSJ,AJ,KR] = caputo_diff(Tstep,GT,GTnow,0,RR, REq, R, U, afd); % Alternative method, not very stable with rebound. Keep for ref.

            S0 = 2*trapz(KR,SJ)/Ze ...
                + (-5+REq^4+4*REq)/(2*Ca);
           
            % (End FDKV model initial stress computation)
            
            display(S0);

        end
        
    elseif zzzen == 1 % Revised Zener model - May. 2021
        % For Maxwell element, we have exact solution at t=0. Just set up storage structure for later steps. 
        
        % Similar to old Zener model, get stretch at end of expansion:
        lam0 = 1/REq; % Stretch at r = R, t = 0;
        texp = -IX(1,1); % Effective duration of expansion
        tau = mu/G1; % Inherent time scale in Maxwell element;
        
        % We need to perform spatial integration to get stress integral in
        % Maxwell element. 
        
        % Borrow spatial mesh from FD model
        RR = 1 + RMesh; % Note: mesh stored as (r0-1) in ref. config.
        RR = RR*REq; % Multiply everything by actual bubble radius
        nr = length(RR);
        
        [SR,RT] = maxwell_stress(RR,lam0,texp,tau); % Ignore modulus for now
        Sneq0 = (1/De)*(1/Re)*trapz(RT,SR);
        
        % Then add ground-state elastic part to get total:
        S0 = -1/(2*Ca)*(5-REq^4-4*REq) + Sneq0;
        
        % Take care of storage structure: (Imitate FD for now. Perhaps not
        % optimized ...)
        nchunk = 5000;
        BT = zeros(nchunk,nr); % 'Beta' at each node. Only use most recent and only for t>=0, but store everything.
        
        % 6/5/21: Note that 'Beta' itself includes exponent of a positive
        % number -- may blow up when t>>mu/G1. To resolve this issue, we
        % will store "Beta*exp(-t*G1/mu)" in BT. For t=0, the exponent is 1. 
        
        % Fill the first row:
        BT(1,:) = (-1/12)*(1/De)*(1/Re)*(SR').*(RT');
        
    elseif fdmax == 1
        % Set up the approximation parameters here:
        [ReK,DeK] = fdtrans(afd,Ze);
        
        % Get size:
        nK = length(DeK);
        
        % Now, evaluate stress in Maxwell sense --
        % (Comments provided for zzzen, above. nhzen is Jon's original version, 
        % with minor revisions by Zhiren.)
        
        % Since Maxwell solver was set up in absolute scale, get absolute
        % relaxation for each element:
        tauK = DeK*t0;
        
        lam0 = 1/REq; % Stretch at r = R, t = 0;
        texp = -IX(1,1); % Effective duration of expansion
        
        RR = 1 + RMesh; % Note: mesh stored as (r0-1) in ref. config.
        RR = RR*REq; % Multiply everything by actual bubble radius
        nr = length(RR);
        
        SNQ = zeros(nK,1); % Set up a matrix to store non-eq. stress for each element
        dSNQ = zeros(nK,1); % Initialize here, even though not necessary ...
        
        % Set up storage earlier, since we have multiple Maxwell elements: 
        % (To speed up computation, let's only store recent value, rather
        % than entire history ... i.e., no more "nchunk")
        BT = zeros(nK,nr); % 'Beta' at each node. 
        
        for kk = 1:nK % Since we need to evaluate stress numerically, this loop is needed
            [SR,RT] = maxwell_stress(RR,lam0,texp,tauK(kk));
            SNQ(kk) = (1/DeK(kk))*(1/ReK(kk))*trapz(RT,SR);
            
            % Storage:
            BT(kk,:) = (-1/12)*(1/DeK(kk))*(1/ReK(kk))*(SR').*(RT');
        end
        
        % Comment: We are repeating the evaluation of RT. Only stress is
        % decoupled; deformation is same. Can potentially speed up ...
        
        % Bring in ground state and sum:
        S0 = -1/(2*Ca)*(5-REq^4-4*REq) + sum(SNQ);
        
        display(S0); % For debug;
        
    end
    
end

X0 = [R0_star U0_star P0_star S0 Tau0 C0 Tm0];
%X0 = [R0_star U0_star P0_star Tau0 C0 Tm0]; % Try dropping S?

tau_del= [];
tdel=[];
Tdel = [];
Cdel = [];

%************************************************
% March equations in time
%options = odeset('RelTol',1e-10);

trange = linspace(0,tspan_star,1000); % This allows us to force truncation at certain temporal positions, if needed
% Set up temporal checking points:
tplus = 0.3*t0; %3E-5; % absolute time
tcheck = tplus;

% ZZ - use modified ODE solver for FDKV
if fdkv == 1 || zzzen == 1 || fdmax == 1
    % Globalize selected variables:
    %global tout
    
    % New strategy: store 'tout' in .mat file, named
    % 'intermediate_tout_storage.mat'. Store 'tout' inside from solver.
    
    % global rb vb udot    % Maybe udot is already passed? Confirm.
    
    % Instead of globalize, just initialize the variables:
    rb = 0;
    vb = 0;
    udot = 0;
    % tout = [];  % No longer used 
    
    % Initialize storage file: (Works, but slowed down code)
    % save('intermediate_tout_storage.mat','tout');
    
    % Better approach: use OutputFcn in ode23tb
    opts0 = odeset('OutputFcn',@trackstep, 'RelTol',1e-8,'AbsTol',1e-8); % Define option for ODE solver
    
    % Run
    %[t , X] = ode23tb_spit(@bubble, trange, X0,opts0);
    [t , X] = ode23tb(@bubble, trange, X0, opts0);
else
    opts1 = odeset('RelTol',1e-8,'AbsTol',1e-8); % Define option for ODE solver
    [t , X] = ode23tb(@bubble, trange, X0, opts1);
    %[t , X] = ode23tb_spit(@bubble, trange, X0); % See if this messes things up
end

% Clear the global variables: (2/27 - Try doing this in run file)
% clear tout rb vb udot

R = X(:,1); % Bubble wall Radius
U = X(:,2); % Bubble wall velocity
P = X(:,3); % Internal pressure
S = X(:,4); % Stress integral
Tau = X(:,5:(NT+4)); % Variable relating to internal temp   - Alt: X(:,4:(NT+3)); 
C =  X(:,(NT+5):(2*NT+4)); % Vapor concentration in the bubble  - Alt: X(:,(NT+4):(2*NT+3)); 
Tm = X(:, (2*NT+5):end ); % Temperature variation in the medium - Alt: X(:, (2*NT+4):end ); 
T = (A_star -1 + sqrt(1+2*Tau*A_star)) / A_star; % Temp in bubble

% ******************************
% Transform variables back into their dimensional form
if (Dim == 1)
    R = R*R0;
    t = t*t0;
    T = T*T_inf;
    P = P*P_inf;
    S = []; %S*P_inf;
    U = U*(R0/t0);
    tdel= tdel*t0;
    Tdel = Tdel*T_inf;
end
%***********************


%*************************************************************************
% Nested function; ODE Solver calls to march governing equations in time
% This function has acess to all parameters above

    function dxdt = bubble(t,x)
        
        % ZZ - drop stress integral from update
        
        % Break x vector into indv. values
        R = x(1); % Bubble wall Radius
        U = x(2); % Bubble wall velocity
        P = x(3); % Internal pressure
        %Z1 = x(4); % Stress integral
        S = x(4);
        Tau = x(5:(NT+4)); % Alt - x(4:(NT+3));
        C = x((NT+5):(2*NT+4)); % Alt - x((NT+4):(2*NT+3)); 
        Tm = x((2*NT+5):end); % Alt - x((2*NT+4):end); 
        
        if (disptime == 1)
            disp(t/tspan_star);
        end
        
        % *********Solves for boundary condition at the wall**************
        %         if (Tmgrad == 1)
        %             if t/tspan_star> 0.001
        %                 %Might need to tune 0.001 for convergence:
        %                 guess= -.001+tau_del(end);
        %                 prelim  = fzero(@Boundary,guess);
        %             else
        %                 guess = -.0001;
        %                 prelim  = fzero(@Boundary,guess);
        %             end
        %         else
        prelim = 0;
        %        end
        
        %****************************************************************
        % Sets value at boundary conditions
        tau_del = [tau_del prelim];
        Tau(end) = prelim;
        T = TW(Tau);
        Tm(1) = T(end);
        % TW(prelim)
        
        % Calculated variables
        K_star = A_star*T+B_star;
        C(end) =  CW(T(end),P);
        
        Rmix = C*Rv_star + (1-C)*Ra_star;
        
        % Gets variables that are not directly calculated as outputs
        Tdel = [Tdel T(end)];
        tdel = [tdel t];
        Cdel = [Cdel C(end)];
        
        %Set external pressure
        %         if (Pext_type == 'sn')
        %             Pext =  -Pext_Amp_Freq(1)/P_inf*sin( Pext_Amp_Freq(2)*t*t0) ;
        %             P_ext_prime = -Pext_Amp_Freq(2)*t0*Pext_Amp_Freq(1)/P_inf...
        %                 *cos( Pext_Amp_Freq(2)*t*t0) ;
        %         elseif (Pext_type == 'RC')
        %             Pext = Pext_Amp_Freq(1)/P_inf ;
        %             P_ext_prime = 0;
        %         elseif (Pext_type == 'RG')
        %             Pext = -Pext_Amp_Freq(1)/P_inf ;
        %             P_ext_prime = 0;
        %         elseif (Pext_type == 'ip')
        %             Pext = -Pext_Amp_Freq(1)/P_inf*...
        %                 (1-heaviside(t-Pext_Amp_Freq(2)/t0)) ;
        %             P_ext_prime = 0;
        if strcmp(Pext_type,'IC')
            Pext = 0;
            P_ext_prime = 0;
        elseif strcmp(Pext_type , 'ga')
            Pext = pf(t)/P_inf;
            P_ext_prime = pfdot(t)/P_inf;
        end
        
        % *****************************************
        % Create derivative terms
        
        % Temp. field of the gas inside the bubble
        DTau  = D_Matrix_T_C*Tau;
        DDTau = DD_Matrix_T_C*Tau;
        
        % Concentration of vapor inside the bubble
        DC  = D_Matrix_T_C*C;
        DDC = DD_Matrix_T_C*C;
        
        % Temp. field in the material
        DTm = D_Matrix_Tm*Tm;
        DDTm = DD_Matrix_Tm*Tm;
        
        %***************************************
        % Internal pressure equation
        pdot = 3/R*(Tgrad*chi*(k-1)*DTau(end)/R - k*P*U +...
              + Cgrad*k*P*fom*Rv_star*DC(end)/( T(end)*R* Rmix(end)* (1-C(end)) ) );
        % *****************************************
        
        %***************************************
        % Temperature of the gas inside the bubble
        U_vel = (chi/R*(k-1).*DTau-yk*R*pdot/3)/(k*P);
        first_term = (DDTau.*chi./R^2+pdot).*(K_star.*T/P*(k-1)/k);
        second_term = -DTau.*(U_vel-yk*U)./R;
        
        Tau_prime = first_term+second_term;
        Tau_prime(end) = 0;
        Tau_prime = Tau_prime*Tgrad;
        % *****************************************
        
        %***************************************
        % Vapor concentration equation
        U_mix = U_vel + fom/R*((Rv_star - Ra_star)./Rmix).*DC;
        one = DDC;
        two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
        three =  (U_mix-U.*yk)/R.*DC;
        
        C_prime = fom/R^2*(one - two) - three;
        C_prime(end) = 0;
        C_prime = C_prime*Cgrad;
        %*****************************************
        
        %***************************************
        % Material temperature equations
        first_term = (1+xk).^2./(L*R).*(U./yk2.^2.*(1-yk2.^3)/2+foh/R.*((xk+1)/(2*L)-1./yk2)).* DTm;
        second_term = foh/R^2.*(xk+1).^4/L^2.*DDTm/4;
        third_term =  3*Br./yk2.^6.*(4/(3*Ca).*(1-1/R^3)+4.*U/(Re.*R)).*U./R;
        Tm_prime = first_term+second_term+third_term;
        Tm_prime(end) = 0; % Sets boundary condition on temp
        Tm_prime(1) = 0; % Previously calculated;
        Tm_prime = Tm_prime*Tmgrad; %Tmgrad makes this quantity zero
        %*****************************************
        
        % ZZ: For debug purpose, display time:
        % t_abs = t*t0;
        % if t_abs > tcheck
        %     display(t_abs);
        %     tcheck = tcheck + tplus; % advance check point
        % end
%         display(t_abs); % For major debug -- e.g., when code is stuck.
%         
        %***************************************
        % Elastic stress in the material
         if linkv == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -4/(3*Ca)*(1 - 1/Rst^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R/Rst^3 + 4/Re*U^2/R^2;
            else
                S = -4/(3*Ca)*(1 - 1/R^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R^4 + 4/Re*U^2/R^2;
            end
        elseif neoHook == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re*U/R ;
                Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re*U^2/R^2;
            else
                S = -(5 -4/R - 1/R^4)/(2*Ca) - 4/Re*U/R;
                Sdot =  -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re*U^2/R^2;
            end
        elseif Yeoh == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) ...
                    + 4/CaY*(177/40+Rst.^-8/8+Rst.^-5/5-3*Rst.^-4/4+Rst.^-2-3./Rst-2*Rst) ...
                    - 4/Re*U/R;
                Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca...
                    + 4/CaY*U/R*(-Rst.^-8-Rst.^-5+3*Rst.^-4-2*Rst.^-2+3./Rst-2*Rst)...
                    + 4/Re*U^2/R^2;
            else
                S = -(5 -4/R - 1/R^4)/(2*Ca)...
                    + 4/CaY*(177/40+R.^-8/8+R.^-5/5-3*R.^-4/4+R.^-2-3./R-2*R) ...
                    - 4/Re*U/R;
                Sdot =  -2*U*(1/R^2 + 1/R^5)/Ca...
                    +4/CaY*(-R.^-9-R.^-6+3*R.^-5-2*R.^-3+3*R^-2-2)...
                + 4/Re*U^2/R^2;
            end
        elseif nhkv_pld == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                %if abs(U)>0
                %    pause(0.01)
                %end
                %S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re*sign(U)*(2^alpha+1)/(3*alpha)*(abs(U)/R)^alpha;
                S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re*sign(U)*(2^alpha+1)/(3*alpha)*(abs(U)/R)^(alpha);
                %Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re/3*(2^alpha+1)*(abs(U)/R)^(alpha+1);
                Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re/3*(2^alpha+1)*sign(U)*(abs(U)/R)^(alpha);
            else
                S = -(5 -4/R - 1/R^4)/(2*Ca) - 4/Re*sign(U)*(2^alpha+1)/(3*alpha)*(abs(U)/R)^alpha;
                Sdot =  -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re/3*(2^alpha+1)*(abs(U)/R)^(alpha+1);
            end
        elseif sls == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                Sdot = -S/De - 4*(1-1/Rst^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            else
                Sdot = -S/De - 4*(1-1/R^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            end
        elseif nhzen == 1
            if  strcmp(Pext_type, 'IC')
                Rst = R/REq;
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                if isinf(Sdot)
                    Rst=Rst+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                end
            else
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                if isinf(Sdot)||isnan(Sdot)
                    R = R+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                end
            end
            
        elseif zzzen == 1
            % Revised Zener formulation -- ZZ, May. 2021
            % Again, imitate FD code structure for now
            
            tspit = [0,trackstep([],[],'get')]; % Add t=0 at beginning, since it does not get added
            tnow = tspit(end);
            
            nstep = length(tspit);
            
            % Deal with t = 0 separately:
            if nstep == 1
                S = S0;
                Sdot = -Sneq0/De;
            else
                tprev = tspit(end-1);
                tdiff = tnow - tprev;
                
                % Get current deformation information:
                RK = ((R^3-REq^3) + RR.^3).^(1/3);
                VK = U*(R^2)*(RK).^(-2);
                Bprev = BT(nstep-1,:);
                % Bnow = Bprev + (VK./RK)'*De*(exp(tnow/De)-exp(tprev/De));  % This form produced singularity. Bring exp(-tnow/De) in from outside and store Beta*exp(-t/De).
                Bnow = Bprev*exp(-tdiff/De) + (VK./RK)'*De*(1-exp(-tdiff/De));
                
                % Save the values, even if we end up replacing it:
                BT(nstep,:) = Bnow;
                
                % Take intermediate steps, just to make life easier
                %SRR = -4*exp(-tnow/De)*Bnow/(De*Re);              % Radial stress - Abandoned form due to beta singularity
                SRR = -4*Bnow/(De*Re);
                DS = -4*(VK./RK)'/(De*Re)-SRR/De;    % Rate of radial stress
                SI1 = 3*SRR./(RK'); % integrand for S
                SI2 = 3*DS./(RK') - 3*SRR.*(VK./(RK.^2))'; % integrand for Sdot
                
                % Integrate:
                Sneq = trapz(RK,SI1);
                Sdneq = trapz(RK,SI2);
                    
                if strcmp(Pext_type, 'IC')
                    Rst = R/REq;
                else
                    Rst = R;    
                end 

                S = Sneq + (-5+Rst^(-4)+4/Rst)/(2*Ca);
                Sdot = Sdneq - 2*(U/R)*(Rst^(-4)+1/Rst)/Ca;
            end
            
        elseif fdmax == 1   % Maxwell approximation of FD, Dec. 2021
            % Procedure similar to zzzener, just more elements using same
            % update algorithm.
            
            tspit = [0,trackstep([],[],'get')]; % Add t=0 at beginning, since it does not get added
            tnow = tspit(end);
            
            nstep = length(tspit);
            
            % Deal with t = 0 separately:
            if nstep == 1
                S = S0;
                dSNQ = SNQ./(DeK'); 
                Sdot = sum(dSNQ);
            else 
                tprev = tspit(end-1);
                tdiff = tnow - tprev;
                
                % Get current deformation info:
                RK = ((R^3-REq^3) + RR.^3).^(1/3);
                VK = U*(R^2)*(RK).^(-2);
                
                % (Conveniently, BT is already current)
                
                % Note: ReK & DeK are row vectors; RK and VK are column.
                % Diagonalize vector to multiple each row by assigned value
                Bnow = diag(exp(-tdiff./DeK))*BT + (DeK.*(1-exp(-tdiff/De)))'*(VK./RK)';
                BT = Bnow; % Update
                
                % Take intermediate steps
                SRR = -4*diag(1./(DeK.*ReK))*Bnow;  % [nK x nr]
                DS = -diag(1./DeK)*SRR - 4*(1./(DeK.*ReK))'*(VK./RK)'; % [nK x nr]
                SI1 = 3*SRR*diag(1./RK);
                SI2 = 3*DS*diag(1./RK) - 3*SRR*diag(VK./(RK.^2));
                
                % Integrate (Loop to get value for each element)
                for kk = 1:nK
                    SNQ(kk) = trapz(RK,SI1(kk,:));
                    dSNQ(kk) = trapz(RK,SI2(kk,:));
                end
                
                if strcmp(Pext_type, 'IC')
                    Rst = R/REq;
                else
                    Rst = R;    
                end 
                
                S = sum(SNQ) + (-5+Rst^(-4)+4/Rst)/(2*Ca);
                Sdot = sum(dSNQ) - 2*(U/R)*(Rst^(-4)+1/Rst)/Ca;
            end
                

            
        elseif fdkv == 1    % FDKV - ZZ, Feb. 2021
            
            % Update history storage
            
            % Load latest version of tout from storage:
            % load('intermediate_tout_storage.mat');
            
            %[tnow,loc] = max(tout); % Get current time from ODE solver
         
            % Faster method, use option in ODE solver and store data:
            tspit = [0, trackstep([],[],'get')]; % Add t = 0 at beginning - also ensure tspit is not empty
            
            % With new structure, we don't even need to remove zeros!
            tnow = tspit(end);

 
            if tnow > Tstep(end); % Compare with stored steps
                %Tstep = [TX;(tout(1:loc))']; % Note: tout(1) = 0   % Only a vector, not big matrix. I will allow append for convenience ...
                Tstep = [TX;tspit'];
                nstep = length(Tstep);
            
                abk = udot; % Most recent solution of bubble acceleration
                %Astore = [Astore; zeros(1,nr)]; % Old structure, no longer needed due to pre-allocation.
                
                % Instead of loop, do most operation with matrix: (GT is
                % 3D, so leave inside loop ...)
                
                % We have RK, VK already from previous step. (Not RRBJ)
                
                %RK = ((rb^3-REq^3) + RR.^3).^(1/3);
                %VK = vb*(rb^2)*(RK).^(-2);
                RRBJ = rb./RK;
                Astore(:,kk) = abk*(RRBJ.^2) + 2*(VK.^2).*RRBJ./RK.*(1-RRBJ.^3); 
                
                if afd > 0 % Just to allow pure elasticity option to skip
                    for jj = 1:nr
                        for ii = 1:8
                            GT(nstep,jj,ii) = GTnow(ii,jj);
                        end

                        %rj0 = RR(jj);   % Undeformed radius at node
                        %rjk = (rb^3+rj0^3-REq^3)^(1/3); % Current radius at node
                        %vjk = vb*(rb/rjk)^2;  % Current velocity at node
                        %rrbj = rb/rjk; % Current radii ratio between bubble wall & node
                        %Astore(nstep,jj) = abk*rrbj^2+2*(vb^2)*rrbj/rjk*(1-rrbj^3);  % Note that Astore is used only for "exact derivative method" of FD eval.

                        % Note: rb & vb were snatched before R & U were updated 
                    end % j - spatial nodes
                end
            end % Time step advancement check
            
            % Evaluate GTnow (We might be stepping back & forth locally, so evaluating each time is safer.)
            
            % Do matrix operation instead of loop:
            RK = ((R^3-REq^3) + RR.^3).^(1/3);
            VK = U*(R^2)*(RK).^(-2); 
            
            for jj = 1:nr
                %rj0 = RR(jj);   % Undeformed radius at node
                %rjk = (R^3+rj0^3-REq^3)^(1/3);
                %vjk = U*(R/rjk)^2;
                rjk = RK(jj);
                vjk = VK(jj);

                GTnow(:,jj) = evalgt(rjk,vjk);
            end % j - spatial nodes
            
            % Evaluate stress integral:
            
            if strcmp(Pext_type, 'IC')
                Rst = R/REq;
            else
                Rst = R;    % ZZ: Double check if this is correct. Did not check carefully since we are interested in Flynn collapse.
            end 

            if afd == 0
                S = -(1/(2*Ca)+1/(2*Ze))*(5-Rst^(-4)-4/Rst);
                Sdot = -2*(U/R)*(Rst^(-4)+1/Rst)*((1/Ca)+(1/Ze));
                Sdd = 0;
            else
                [SJ,dSJ,AJ,KR] = caputo_exact(Tstep,GT,GTnow,t, RR, REq, R, U, Astore, afd); 
                S = 2*trapz(KR,SJ)/Ze ...
                    + (-5+Rst^(-4)+4/Rst)/(2*Ca);
                Sdot = 2*trapz(KR,dSJ)/Ze ...
                    - 2*(U/R)*(Rst^(-4)+1/Rst)/Ca;
                Sdd = 2*trapz(KR,AJ)/Ze; % Note - no contribution to elastic spring
                
                % Sdd is SdotA. Just pass it on for now.
            end
        end
        %****************************************************
        
        % Equations of motion
        rdot = U;
        
        if (Tgrad == 0)
            pdot = -3*k*U/R*P;
        end
        
        Pv = (Pvsat(T(end)*T_inf)/P_inf);
        if comp == 0
            %Rayleigh-Plesset equation
            udot = (P + abs(1-Cgrad)*Pv  - 1 - Pext + S - 1/(We*R) -1.5*U^2)/R;
        else
            % Keller-Miksis equation
            if linkv==1 || neoHook==1 || Yeoh==1
                SdotA = 4/Re;
            elseif sls==1 || nhzen==1 || fdkv==1 || zzzen==1 || fdmax==1  % ZZ - Not exactly true for FDKV, but I also don't think this matters ...
                SdotA = 0;
            elseif nhkv_pld==1
                %SdotA = 4/Re*(2^alpha+1)/3*(abs(U)/R)^(alpha-1);
                SdotA = 4/Re/3*(2^alpha+1)*sign(U)*(abs(U)/R)^(alpha)*R^2/U^2;
                if isnan(SdotA)
                    SdotA=4/Re;
                end
            end
            
            if fdkv == 1
                RHS = (1+U/C_star)...
                    *(P  + abs(1-Cgrad)*Pv -1/(We*R) + S - 1 - Pext)  ...
                    + R/C_star*(pdot+ U/(We*R^2) + Sdot -P_ext_prime );
                LHS = (3/2)*(1-U/(3*C_star))*U^2;
                denom = (1-U/C_star)*R - (R/C_star)*Sdd;
                
                udot = (RHS - LHS)/denom;
                
            else  % Original expression
                 udot = ((1+U/C_star)...
                    *(P  + abs(1-Cgrad)*Pv -1/(We*R) + S - 1 - Pext)  ...
                    + R/C_star*(pdot+ U/(We*R^2) + Sdot - P_ext_prime ) ...
                    - 1.5*(1-U/(3*C_star))*U^2)/((1-U/C_star)*R); % +JdotA/(C_star));
            end
        end
        % ****************************************
        Sdot = Sdot - SdotA*udot/R;     % We don't really use advancement of S in next step though ...
        dxdt = [rdot; udot; pdot; Sdot; Tau_prime; C_prime; Tm_prime];
        
        %dxdt = [rdot; udot; pdot; Tau_prime; C_prime; Tm_prime]; % Try dropping S
        
        % ZZ - Store current R & U globally for FDKV before advancing:
        if fdkv == 1
            rb = R;
            vb = U;
        end
        
    end
%*************************************************************************

% Other nested functions used to carry out repetetive calculations
% throughout the code
    function Tw = TW(Tauw)
        %calculates the temperature at the bubble wall as a fuction of \tau
        Tw = (A_star -1 + sqrt(1+2*Tauw*A_star)) / A_star;
    end

    function Cw = CW(Tw,P)
        % Calculates the concentration at the bubble wall
        
        %Function of P and temp at the wall
        theta = Rv_star/Ra_star*(P./(Pvsat(Tw*T_inf)/P_inf) -1);
        Cw = 1./(1+theta);
    end

    function Tauw = Boundary(prelim)
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

% Gaussian pressure functions
% acoustic pressure
    function p = pf(t)
        if t<(dt_star-5*tw_star) || t>(dt_star+5*tw_star)
            p=0;
        else
            p = -Pext_Amp_Freq(1)*exp(-(t-dt_star).^2/tw_star^2);
        end
    end

% time derivative of acoustic pressure
    function pdot = pfdot(t)
        if t<(dt_star-5*tw_star) || t>(dt_star+5*tw_star)
            pdot=0;
        else
            pdot = 2*(t-dt_star)/tw_star^2*Pext_Amp_Freq(1).*exp(-(t-dt_star).^2/tw_star^2);
        end
        
    end
end

% ---------------------------------------------------------------
% ZHIREN'S FUNCTIONS FOR FD EVALUATION:
% ---------------------------------------------------------------
% (A) FD-related term evaluator
% 8/16/20: We have simpler functions now. Woot.
function GG = evalgt(rjk,vjk)
    %lam = rjk/rj0; % No longer used. We cancelled out the rj0 values.
    
    % Note: RR is already scaled to non-dimensional
    
    g1 = 1/rjk^2;
    g2 = rjk^4;
    g3 = vjk*rjk^(-3);
    g4 = vjk*rjk^(3);
    g5 = 3*vjk^2/rjk^4; 
    g5a = rjk^(-3); % This term is associated with acceleration 
    g6 = 3*vjk^2*rjk^2;
    g6a = rjk^3; % This term is associated with acceleration 
    
    GG = [g1;g2;g3;g4;g5;g5a;g6;g6a];
    
end % Sub-function (A): evalgt

% ---------------------------------------------------------------
% (B1) Caputo with Finite Difference Approach
% [Have not been maintained since summer 2020. Copied here just for reference.]
function [SJ,dSJ,AJ,KR] = caputo_diff(Tstep,GT,GTnow,t, RR, REq, R, U, afd)

% Though no longer used, revise input structure similar to option (B2),
% transfer function to local instead of nested.
% Include global variables inputs: afd, RR, R, REq, U
% Get nr from size of RR
% Astore not used here

    SJ = zeros(1,nr);
    dSJ = zeros(1,nr);
    AJ = zeros(1,nr); % Not really used for this method, but return it blank
    KR = zeros(1,nr); % Current position - remember to integrate over this, not reference!
    
    nstep = length(Tstep);
    
    for jj = 1:nr
        DGJ = zeros(4,1); % Only need 4 for this method
        
        for kk = 2:nstep % Reminder: Tstep(2)-Tstep(1) is step 1
            t1 = Tstep(kk-1);
            t2 = Tstep(kk);
            dt = t2 - t1;
           
            cvdt = dt^(-afd)*(((t-t1)/dt)^(1-afd)-((t-t2)/dt)^(1-afd));
            
            for ii = 1:4
                dgi = GT(kk,jj,ii)-GT(kk-1,jj,ii);
                DGJ(ii) = DGJ(ii) + dgi*cvdt;
            end % i - FD-related terms
            
        end % k - time step
        
        % Finally, treat current step
        dt0 = t - Tstep(end);
        cvdt = dt0^(-afd);
        
        for ii = 1:4
            dgi = GTnow(ii,jj)-GT(nstep,jj,ii);     % Reminder: Don't use GT(end ...)
            DGJ(ii) = DGJ(ii) + dgi*cvdt;
        end % i - FD-related terms
        
        % Apply gamma function altogether
        DGJ = DGJ/gamma(2-afd);
        
        % Retrieve node kinematic info and assemble integrand
        rj0 = RR(jj);                   % Undeformed radius at node
        rjk = (R^3+rj0^3-REq^3)^(1/3);  % Current radius at node
        vjk = U*(R/rjk)^2;            % Current velocity at node
        
        % lam = rjk/rj0; % No longer used
        
        SJ(jj) = rjk*DGJ(1) - DGJ(2)/(rjk^5);
        dSJ(jj) = vjk*(DGJ(1) + 5*DGJ(2)/(rjk^6)) - 2*rjk*DGJ(3) - 4*DGJ(4)/(rjk^5);
        % AJ not evaluated for this method.
        KR(jj) = rjk;
    end % j - spatial node

end % Sub-function (B1): caputo_diff

% Note:
% This method does not perform very well. Reason being that acceleration is
% not properly considered. Keep here just for reference. Use (B2) instead for FD evaluation.

% ---------------------------------------------------------------
% (B2) Caputo with "Exact" First Derivative
function [SJ,dSJ,AJ,KR] = caputo_exact(Tstep,GT,GTnow,t, RR, REq, R, U, Astore, afd)
    
% Updated Notes: 2/25/2021
% Include global variables inputs: afd, Astore, RR, R, REq, U
% Get nr from size of RR
    
    dbstop if error % For debug
    
    nr = length(RR);
    
    SJ = zeros(1,nr);
    dSJ = zeros(1,nr);
    AJ = zeros(1,nr); 
    KR = zeros(1,nr); % Current position - remember to integrate over this, not reference!
    
    nstep = length(Tstep);
    
    %bmark = 0;
    
    % Take care of (most of) CV outside:
    CV0 = (t-Tstep).^(1-afd);
    DCV =  diff(CV0);
    
    for jj = 1:nr
        DGJ = zeros(6,1); % Rates of Q1-Q4 + Two current acceleration terms
        
        % 2/15/21: ODE solver sometimes step back. Need to drop those cases.
        
        for kk = 2:nstep % Reminder: Tstep(2)-Tstep(1) is step 1
            %t1 = Tstep(kk-1);
            t2 = Tstep(kk);
            
            % Check if t2<t indeed
            if t2<t
                %cv = (t-t1)^(1-afd)-(t-t2)^(1-afd);
                cv = -DCV(kk-1);

                % For this method, we need to update the terms separately ...
                DGJ(1) = DGJ(1) - cv*(2*GT(kk,jj,3));
                DGJ(2) = DGJ(2) + cv*(4*GT(kk,jj,4));
                DGJ(3) = DGJ(3) + cv*(-GT(kk,jj,5) + GT(kk,jj,6)*Astore(jj,kk)); % Astore flipped - 2/25/2021
                DGJ(4) = DGJ(4) + cv*(GT(kk,jj,7) + GT(kk,jj,8)*Astore(jj,kk)); % Astore flipped - 2/25/2021
                
                % Save t2 for final step:
                tfin = t2;
            %else
                %bmark = 1;
            end
        end % k - time step
        
        % Finally, treat current step
        %dt0 = t - Tstep(end);
        dt0 = t - tfin; % Use this definition instead
        % Sometimes
        
        cv = dt0^(1-afd);
        
        DGJ(1) = DGJ(1) - cv*(2*GTnow(3,jj));
        DGJ(2) = DGJ(2) + cv*(4*GTnow(4,jj));
        DGJ(3) = DGJ(3) + cv*(-GTnow(5,jj));
        DGJ(4) = DGJ(4) + cv*(GTnow(7,jj)); 
        
        DGJ(5) = cv*GTnow(6,jj); % No previous terms ...
        DGJ(6) = cv*GTnow(8,jj); % No previous terms ...
        
        % Apply gamma function altogether
        DGJ = DGJ/gamma(2-afd);
        
        % Retrieve node kinematic info and assemble integrand
        rj0 = RR(jj);                   % Undeformed radius at node
        rjk = (R^3+rj0^3-REq^3)^(1/3);  % Current radius at node
        vjk = U*(R/rjk)^2;            % Current velocity at node
        
        % lam = rjk/rj0; % No longer used 
        
        SJ(jj) = rjk*DGJ(1) - DGJ(2)/(rjk^5); 
        dSJ(jj) = vjk*(DGJ(1) + 5*DGJ(2)/(rjk^6)) - 2*rjk*DGJ(3) - 4*DGJ(4)/(rjk^5);
        AJ(jj) = -2*rjk*DGJ(5) - 4*DGJ(6)/(rjk^5);
        KR(jj) = rjk;
        
        % 2/15/21 - Finding complex output. Add a break for debug.
        if imag(SJ(jj))~=0
            error('Pause. We are getting complex values.')
        end
        
    end % j - spatial node

%     if bmark > 0
%         display('Solver stepped back'); % Display message, mainly for debug reference.
%     end
    
end % Sub-function (B2): caputo_exact

% Note:
% This method is much more stable. Use this instead of (B1).

% ---------------------------------------------------------------
% (C) Function to track succesful steps in ODE solver
% (Got this from Jan @ MATLAB forum on 2/27/2021)
% (https://www.mathworks.com/matlabcentral/answers/757774-passing-time-step-from-ode-solver-ode23tb-to-ode-function?s_tid=srchtitle)

function status = trackstep(t,y,flag)

persistent tlocal tcount

switch char(flag)
    case ''
        tcount = tcount + 1;
        if tcount > numel(tlocal)
            tlocal(tcount + 1000) = 0; % increase size
        end
        tlocal(tcount) = t;
        status = 0;
        
    case {'init','done'}
        tlocal = zeros(1,1000);
        tcount = 0;     % I tried setting this to 1, seems to upset MATLAB
        status = 0;
        
    case 'get' 
        status = tlocal(1:tcount);
end

end

% Additional reference on setting up own function for ODE solver option:
% https://www.mathworks.com/help/matlab/ref/odeset.html

% ---------------------------------------------------------------
% (D) Initial stress integral processor for Maxwell element in Zener model
function [SR,RT] = maxwell_stress(RR,lam0,texp,tau)

    nr = length(RR);
    rb0 = RR(1);
    rbt = lam0*rb0;
    
    SR = zeros(nr,1);
    RT = zeros(nr,1);
    
    for kk = 1:nr
        rk0 = RR(kk);
        rkt = (rk0^3 + rbt^3 - rb0^3)^(1/3);
        klam = rkt/rk0;
        [~,XM] = maxwell_expand(klam,texp,tau); % Get elastic stretch here
        % (We can get time output & plot to verify if we needed for debug)
        elam = XM(end);
        
        RT(kk) = rkt;
        SR(kk) = -12*log(elam)/rkt;
    end
end

% ---------------------------------------------------------------
% (E) Translate FD to Maxwell approximation
function [ReK,DeK] = fdtrans(afd,Ze)
    % Generate Reynold's & Deborah number for each element
    
    % (1) Time scale is constant in dimensionless space:
    % (We can coarsen for debug purpose)
    zmin = -8;
    zmax = 6;
    gz = 1; %1/8;
    
    DeK = 10.^(zmin:gz:zmax); % This includes "z0", actually
    
    % (2) Get corresponding viscosity:
    scaler = (Ze*gamma(afd)*gamma(1-afd)*afd)/((10^gz)^afd-1);
    corr = 1; %((1+10^gz)/2)^(afd);    % Empirical correction
    
    ReK = corr*scaler*DeK.^(afd-1);
end
% ---------------------------------------------------------------




