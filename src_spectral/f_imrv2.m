function varargout =  f_imrv2(varargin)
% IMR V2
% TODO WRITE INTRO TO THE CODE

% Original Authors: Matthew Warnez & Carlos Barajas
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
%   The various options can be found below.  Code may also run without any
%   inputs.
par = f_call_params(varargin{:});
params = cell2mat(par);
% numerical settings 
polytropic      = params(1); cold = params(2); rayleighplesset = params(3);
enthalpy        = params(4); gil = params(5); neoHook = params(6);
voigt           = params(7); linelas = params(8); liner = params(9);
oldb            = params(10); ptt = params(11); gies = params(12);
% output options
dimensionalout  = params(13); progdisplay = params(14); detail = params(15);
plotresult      = params(16); radiusonly = params(17); vitalsreport = params(18);
displayonly     = params(19); technical = params(20);
% solver options
method          = params(21);
spectral        = params(22);
divisions       = params(23);
% numerical parameters 
Nv              = params(24); Nt = params(25); Mt = params(26);
Lv              = params(27); Lt = params(28);
% physical parameters%
% acoustic parameters
Cstar           = params(29); GAMa = params(30); kappa = params(31);
nstate          = params(32);
% dimensionless waveform parameters
tfin            = params(33); om = params(34); ee = params(35); 
tw              = params(36); dt = params(37); mn = params(38); 
% dimensionless viscoelastic
We              = params(39); Re8 = params(40); DRe = params(41); 
v_a             = params(42); v_nc = params(43); Ca = params(44);
LAM             = params(45); De = params(46); JdotA = params(47); 
v_lambda_star   = params(48); 
% dimensionless thermal 
Foh             = params(49); Br = params(50); alpha = params(51); 
chi             = params(52); iota = params(53);
% dimensionaless mass transfer 
Fom             = params(54); C0 = params(55); Rv_star = params(56);
Ra_star         = params(57); L_heat_star = params(58); mv0 = params(59); 
ma0             = params(60); 
% dimensionless initial conditions
Rzero           = params(61); Uzero = params(62); pzero = params(63);
P8              = params(64); T8 = params(65); Pv_star = params(66);
p0star          = pzero;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical setup and precomputations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collocation point construction
y = cos(pi*(0:Nt)'/(2*Nt));
xi = cos(pi*(0:Mt)'/Mt);
ze = cos(pi*(1:Nv)'/Nv);
% collocation matrix construction
[gA,gAI,~,~,gAPd,gAPdd] = dcdmtxe(Nt);
[mA,~,~,~,mAPd,mAPdd] = dcdmtx(Mt);
Q = [gA(2:end,:) zeros(Nt,Mt+1);
    zeros(Mt,Nt+1) mA(2:end,:);
    2*(0:Nt).^2 iota*(0:Mt).^2];
Q = sparse(Q);
[sCA,sCI,sCAd,~,~,~] = dcdmtx(Nv);
sCA = sCA(2:end,2:end) - 1;
sCI = sCI(2:end,2:end);
sCAd = sCAd(2:end,2:end);
% precomputations
LDR = LAM*De/Re8;
sam = pzero - We + GAMa;
no = (nstate-1)/nstate;
kapover = (kappa-1)/kappa;
yT = 2*Lt./(1+xi) - Lt + 1;
yV = 2*Lv./(1-ze) - Lv + 1;
nn = ((-1).^(0:Nt).*(0:Nt).^2)';
nn = sparse(nn);
Udot = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% precomputations for viscous dissipation
zT = 1 - 2./(1 + (yT - 1)/Lv);
ZZT = cos(acos(zT)*(1:Nv)) - 1;

% compute or load integral coefficients
switch Lv
    case 3, load('eNstore3','eNN'); cdd = eNN(1:Nv)';
    case 4, load('eNstore4','eNN4'); cdd = eNN4(1:Nv)';
    case 5, load('eNstore5','eNN5'); cdd = eNN5(1:Nv)';
    case 6, load('eNstore6','eNN6'); cdd = eNN6(1:Nv)';
    case 10, load('eNstore10','eNN10'); cdd = eNN10(1:Nv)';
    case 20, load('eNstore20','eNN20'); cdd = eNN20(1:Nv)';
    case 0.1, load('eNstorep1','eNNp1'); cdd = eNNp1(1:Nv)';
    case 0.01, load('eNstorep01','eNNp01'); cdd = eNNp01(1:Nv)';
    case 0.3, load('eNstorep3','eNNp3'); cdd = eNNp3(1:Nv)';
    case 0.5, load('eNstorep5','eNNp5'); cdd = eNNp5(1:Nv)';
    otherwise, cdd = preStressInt(Lv,Nv);
end

% index management
if liner == 0
    zeNO = 1; 
else 
    zeNO = 0; 
end
if spectral == 0, Nv = 1; end
if polytropic == 1, Nt = -1; Mt = -1; qdot = []; end
if cold == 1, Mt = -1; end
ia = 4:(4+Nt);
ib = (5+Nt):(5+Nt+Mt);
ic = (6+Nt+Mt):(5+Nt+Mt+Nv);
id = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);
ie = (6+Nt+Mt+2*Nv):(5+2*Nt+Mt+2*Nv);

% initial condition assembly
init = [Rzero; Uzero; pzero; % radius, velocity, pressure
    zeros(Nt+1,1); % auxiliary temperature spectrum
    ones(Mt ~= -1); zeros(Mt,1); % medium temperature spectrum
    zeros(2*(Nv - 1)*(spectral == 1) + 2,1); % stress spectrum
    0]; % initial stress integral
    %zeros(Nt+1,1)]; % initial vapor concentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SOLVER CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = linspace(0,tfin,detail);
% tspan = [0 tfin];
stepcount = 0;

if method == 15
    if divisions == 0
        options = odeset();
    else
        options = odeset('MaxStep',tfin/divisions,'RelTol',1e-4);
    end
    [t,X] = ode15s(@SVBDODE,tspan,init,options);
elseif method == 23
    if divisions == 0
        options = odeset();
    else
        options = odeset('MaxStep',tfin/divisions);
    end
    [t,X] = ode23s(@SVBDODE,tspan,init,options);
elseif method == 45
    if divisions == 0
        options = odeset('NonNegative',1);
    else
        options = odeset('NonNegative',1,'MaxStep',tfin/divisions,'RelTol',1e-12);
    end
    [t,X] = ode45(@SVBDODE,tspan,init,options);
else
    if divisions == 0
        options = odeset('NonNegative',1);
    else
        options = odeset('NonNegative',1,'MaxStep',tfin/divisions);
    end
    [t,X] = ode45(@SVBDODE,tspan,init,options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% functions called by solver %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % acoustic pressure
% mn = 3.7;
% function p = pf(t)
%     
%     if t < dt - pi/om
%         p = 0;
%     elseif t > dt + pi/om
%         p = 0;
%     else
%         p = ee*(0.5 + 0.5*cos(om*(t - dt))).^mn;
%     end
%     
% end
% 
% % time derivative of acoustic pressure
% function pdot = pfdot(t)
%     if t < dt - pi/om
%         pdot = 0;
%     elseif t > dt + pi/om
%         pdot = 0;
%     else
%         pdot = -ee*mn*(0.5+0.5*cos(om*(t-dt))).^(mn-1)*0.5*om.*sin(om*(t-dt));
%     end
%     
% end

% acoustic pressure
function p = pf(t)
    p = ee;
end

% time derivative of acoustic pressure
function pdot = pfdot(t)
    pdot = 0;
end

% stress differentiator
function [trr,dtrr,t00,dt00] = stressdiff(c,d)
    if Nv < 650
        trr = sCA*c;
        dtrr = sCAd*c;
        t00 = sCA*d;
        dt00 = sCAd*d;
    else
        [trr,dtrr] = fctdShift(c);
        [t00,dt00] = fctdShift(d);
    end
end

% stress solver
function s = stresssolve(x)
    if Nv < 650
        s = sCI*x;
    else
        s = fctShift(x);
    end
end

% fast Chebyshev transform
function a = fctShift(v)
    v = v(:);
    v = [0; v; flipud(v(1:Nv-1))];
    a = real(fft(v))/Nv;
    a = [a(2:Nv); a(Nv+1)/2];
end

% fast Chebyshev transform and differentiate
function [v,w] = fctdShift(a)
    M = Nv + 1;
    a = a(:)';
    dd = Nv*[0 a(1:Nv-1) a(Nv)*2 fliplr(a(1:Nv-1))];
    v = ifft(dd);
    v = v(2:M)' - sum(a);
    n2b = (0:M-2).^2.*dd(1:Nv);
    cc = imag(ifft([0:M-2 0 2-M:-1].*dd));
    w = zeros(Nv,1);
    w(1:Nv-1) = csc(pi/Nv*(1:M-2)).*cc(2:Nv);
    w(Nv) = sum((-1).^(1:Nv).*n2b)/Nv + 0.5*(-1)^M*Nv*dd(M);
end

%%%%%%%%%%%%%%%%%%%
% solver function %
%%%%%%%%%%%%%%%%%%%
function dXdt = SVBDODE(t,X)
    
    stepcount = stepcount + 1;
    if progdisplay == 1, disp(t/tfin); end
    
    % extract standard inputs
    R = X(1); U = X(2); p = X(3);

    % non-condensible gas pressure and temperature
    if polytropic == 1
        % polytropic approximation
        p = p0star*R^(-3*kappa);
        pdot = -3*kappa*U/R*p;
        pVap = Pv_star;
    else
        % extract auxiliary temperature
        SI = gA*X(ia);
        % auxiliary temperature derivatives
        dSI = gAPd*SI; % first order derivative
        ddSI = gAPdd*SI; % second order derivative
        % temperature and thermal diffusivity fields
        T = (alpha - 1 + sqrt(1+2*alpha*SI))/alpha;
        D = kapover*(alpha*T.^2 + (1-alpha)*T)/p;
        % new pressure and auxiliary temperature derivative
        pdot = 3/R*((kappa-1)*chi/R*dSI(1) - kappa*p*U);
        
        SIdot = pdot*D + chi/R^2*(2*D./y - kapover/p*(dSI - dSI(1)*y)).*dSI ...
            + chi*D/R^2.*ddSI;
        SIdot(end) = pdot*D(end) - chi/R^2*(8*D(end)*sum(nn.*X(ia)) ...
            + kapover/p*dSI(end)^2) + chi*D(end)/R^2.*ddSI(end);
        
        if cold == 1 % cold-liquid approximation
            % solve auxiliary temperature with boundary condition
            qdot = gAI*[0; SIdot(2:end)];
        else
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
            elseif voigt == 1 || neoHook == 1 || linelas == 1
                TLdot = TLdot + 4*Br./yT.^6*(U/R*(1-1/R^3)/Ca + 3/Re8*(U/R)^2);      
            end
            % enforce boundary condition and solve
            TLdot(end) = 0;
            qdot = [ones(1,Nt+1) -(alpha*(T(1)-1)+1)*ones(1,Mt+1); Q]...
                \[0; SIdot(2:end); TLdot(2:end); 0];       
        end
        pVap = (f_pvsat(T(end)*T8)/P8);
    end
    J = 0; JdotX = 0; Z1dot = 0; Z2dot = 0;
    % stress
    if neoHook == 1 % Kelvin-Voigt with neo-Hookean elasticity
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        % compute stress integral
        J = (4/R + 1/R^4 - 5)/(2*Ca) - 4/Re8*U/R;
        JdotX = -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re8*U^2/R^2;
    elseif voigt == 1 % Kelvin-Voigt (Yang-Church)
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        % compute stress integral
        J = -4/(3*Ca)*(1 - 1/R^3) - 4/Re8*U/R;
        JdotX = -4/Ca*U/R^4 + 4/Re8*U^2/R^2;
    elseif linelas == 1 % Linear elastic
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;     
        % compute stress integral
        J = -2/Ca*(1 - 1/R^2) - 4/Re8*U/R;
        JdotX = -4/Ca*U/R^3 + 4/Re8*U^2/R^2;
    elseif spectral == 1 % Giesekus, PTT, or forced spectral         
        % extract stress spectrum
        c = X(ic); d = X(id);
        % inverse Chebyshev transforms and derivatives
        [trr,dtrr,t00,dt00] = stressdiff(c,d);
        % new spectral coefficient derivatives
        exptau = exp(ptt*Re8*De*(trr + 2*t00));            
        Z1dot = stresssolve(-(exptau/De + zeNO*4*U./(yV.^3*R) ...
            + gies*Re8*trr).*trr ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dtrr ...
            - 4./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re8*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*4*LAM/Re8*(U/R)^2./yV.^6);
        Z2dot = stresssolve(-(exptau/De - zeNO*2*U./(yV.^3*R) ...
            + gies*Re8*t00).*t00 ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dt00 ...
            + 2./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re8*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*10*LAM/Re8*(U/R)^2./yV.^6);
        % compute stress integral
        J = 2*sum(cdd.*(c-d));
        JdotX = 2*sum(cdd.*(Z1dot - Z2dot));   
    elseif liner == 1 % linear Maxwell, linear Jeffreys, linear Zener
        % extract
        Z1 = X(ic);
        J = Z1/R^3 - 4*LAM/Re8*U/R;
        % stress integral derivative
        Z1dot = -Z1/De + 4*(LAM-1)/(Re8*De)*R^2*U - 4*(R^3-1)/(3*Ca*De);
        Z2dot = 0;
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re8*U^2/R^2;     
    elseif oldb == 1 % upper-convected Maxwell, OldRoyd-B
        % extract stress sub-integrals
        Z1 = X(ic); Z2 = X(id);
        % compute new derivatives
        Z1dot = -(1/De - 2*U/R)*Z1 + 2*(LAM-1)/(Re8*De)*R^2*U;
        Z2dot = -(1/De + 1*U/R)*Z2 + 2*(LAM-1)/(Re8*De)*R^2*U;
        J = (Z1 + Z2)/R^3 - 4*LAM/Re8*U/R;
        JdotX = (Z1dot+Z2dot)/R^3 - 3*U/R^4*(Z1+Z2) + 4*LAM/Re8*U^2/R^2;
    end
    
%     %***************************************
%     % Vapor concentration equation 
%     
%     U_mix = U_vel + fom/R*((Rv_star - Ra_star)./Rmix).*DC  ;
%     one = DDC;
%     two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
%     three =  (U_mix-U.*yk)/R.*DC; 
%     % *****************************************
%     
%     % *****************************************
%     % C_prime
%     
%     C_prime = fom/R^2*(one - two) - three;
%     C_prime(end) = 0;
%     C_prime = C_prime*Cgrad;
%     %***************************************
    
    % bubble wall acceleration
    % Rayleigh-Plesset        
    if rayleighplesset == 1
        Udot = (p - 1 + pVap - pf(t) - We/R + J - 1.5*U^2)/R;
    % Keller-Miksis in enthalpy
    elseif enthalpy == 1    
        hB = sam/no*(((p - We/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - We/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/Cstar)*(hB - pf(t)) - R/Cstar*pfdot(t) ...
            + R/Cstar*hH*(pdot + We*U/R^2 + JdotX) ...
            - 1.5*(1 - U/(3*Cstar))*U^2)/((1 - U/Cstar)*R + JdotA*hH/Cstar);
    % Gilmore equation
	elseif gil == 1
        hB = sam/no*(((p - We/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - We/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/Cstar)*(hB - pf(t)) - R/Cstar*pfdot(t) ...
            + R/Cstar*hH*(pdot + We*U/R^2 + JdotX) ...
            - 1.5*(1 - U/(3*Cstar))*U^2)/((1 - U/Cstar)*R + JdotA*hH/Cstar);
    % Keller-Miksis in pressure        
    else    
        Udot = ((1+U./Cstar)*(p - 1 + pVap - pf(t) - We./R + J) ...
            + R./Cstar.*(pdot + We.*U./R.^2 + JdotX - pfdot(t)) ...
            - 1.5.*(1-U./(3.*Cstar)).*U.^2)./((1-U./Cstar).*R + JdotA./Cstar);
    end
    % output assembly
    Jdot = JdotX - JdotA*Udot/R;
    dXdt = [U; Udot; pdot; qdot; Z1dot; Z2dot; Jdot];
end

%%%%%%%%%%%%%%%%%%%
% post-processing %
%%%%%%%%%%%%%%%%%%%

% extract result
R = X(:,1); U = X(:,2); p = X(:,3);
a = X(:,ia)'; b = X(:,ib)'; c = X(:,ic)'; d = X(:,id)';
I = X(:,end);
pA = zeros(size(t));
% Udot = R2dot(t,X);
for n = 1:length(t)
    pA(n) = pf(t(n)); 
end 

% transform to real space
if spectral == 1
    trr = sCA*c; t00 = sCA*d;
else
    trr = c; t00 = d;
end
if polytropic == 0
    T = (alpha-1+sqrt(1+2*alpha*gA*a))/alpha;
    if cold == 0
        TL = mA*b; 
    end
else
    T = R.^(-3*kappa);
end

if dimensionalout == 1
    
    % re-dimensionalize problem
    t = t*tc; R = R*Rref; U = U*uc; p = p*p0; pA = pA*p0; I = I*p0; T = T*T8;
    c = c*p0; d = d*p0; Udot = Udot*uc/tc;
    if spectral == 1
        trr = trr*p0; 
        t00 = t00*p0; 
    end
    if polytropic == 0
        if cold == 0
            TL = TL*T8; 
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%
% output and display %
%%%%%%%%%%%%%%%%%%%%%%

% display figure
if plotresult == 1
    if vitalsreport == 0
        if radiusonly == 1
            plot(t,R); 
            hold('on'); 
            grid('on');
            axis([0 t(end) 0 (max(R) + min(R))]);
        else
            subplot(2 + spectral,1,1);
            plot(t,R); 
            grid('on'); 
            ylabel('R');
            axis([0 t(end) 0 (max(R) + min(R))]);
            subplot(2 + spectral,1,2);
            plot(t,I); 
            grid('on');
            ylabel('J');
            if spectral == 1
                subplot(3,1,3);
                plot(t,trr(end,:),t,t00(end,:));
                ylabel('\tau_{rr}|_R, \tau_{\theta\theta}|_R');
            end
            xlabel('t');
        end
    else
        subplot(3 - polytropic,1,1);
        plot(t,R); 
        hold('on'); 
        grid('on'); 
        ylabel('R');
        axis([0 t(end) 0 (max(R) + min(R))]);
        subplot(3 - polytropic,1,2);
        semilogy(t,abs(c(end,:)),t,abs(d(end,:))); 
        ylabel('c_P, d_P');
        axis([0 t(end) 1e-20 1]);
        if polytropic == 0
            if cold == 1, b = zeros(size(a)); end
            subplot(3,1,3);
            semilogy(t,abs(a(end,:)),t,abs(b(end,:))); 
            ylabel('a_N, b_M');
            axis([0 t(end) 1e-20 1]);
        end
        xlabel('t');
    end
    
end

% assemble output
if displayonly == 1
    varargout = {};
else
    % standard outputs
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = p;
    varargout{4} = U;
    varargout{5} = trr;
    varargout{6} = t00;
    varargout{7} = I;
    varargout{8} = T;    
    if polytropic == 0 && cold == 0
        varargout{9} = TL;
    else
        varargout{9} = ((T8 - 1)*dimensionalout + 1)*ones(divisions,1);
    end   
    % technical data
    if technical == 1
        varargout{10} = a; 
        varargout{11} = b;
        varargout{12} = c; 
        varargout{13} = d;
        varargout{14} = [stepcount p0 alpha];
        varargout{15} = pA;
%         varargout{16} = Udot;
    end 
end

% convert run settings to strings
if rayleighplesset == 1
    eqn = 'Rayleigh Plesset equation';
elseif enthalpy == 1
    eqn = 'Keller-Miksis in enthalpy';
else
    eqn = 'Keller-Miksis in pressure';
end
const = 'none';
if voigt == 0
    if neoHook == 1
        if Ca == Inf
            const = 'Newtonian fluid';
        else
            const = 'neo-Hookean Voigt';
        end
    elseif linelas == 1
        if Ca == Inf
            const = 'Newtonian fluid';
        else
            const = 'linear elastic Voigt';
        end
    elseif liner == 1
        if Ca ~= Inf && LAM == 0
            const = 'linear Zener';
        elseif Ca == Inf && LAM == 0
            const = 'linear Maxwell';
        elseif Ca == Inf && LAM ~= 0
            const = 'linear Jeffreys';
        else
            const = 'Kelvin-Voigt series';
        end
    elseif liner == 0
        if ptt == 0 && gies == 0
            if Ca ~= Inf && LAM == 0
                const = 'upper-convective Zener';
            elseif Ca == Inf && LAM == 0
                const = 'upper-convective Maxwell';
            elseif Ca == Inf && LAM ~= 0
                const = 'Oldroyd-B';
            end
        elseif Ca == Inf && LAM == 0
            if ptt == 1 && gies == 0
                const = 'Phan-Thien-Tanner';
            elseif ptt == 0 && gies ~= 0
                const = ['Giesekus(' num2str(gies) ')'];
            end
        end
    end
elseif liner == 0 && gies == 0 && ptt == 0 && LAM == 0, const = 'Yang-Church Voigt';
end

if polytropic == 0
    if cold == 0, therm = 'full';
    else
        therm = 'cold-medium approximation';
    end
else
    therm = 'polytropic approximation';
end

if spectral == 1
    solut = 'spectral method';
else
    solut = 'ODE formulation';
end

% display run settings
disp('--- Game settings ---');
disp(['Radial dynamics: ' eqn]);
disp(['Medium rheology: ' const]);
disp(['Thermal effects: ' therm]);
disp(['Solution method: ' solut]);
disp('--- Dimensionless numbers ---');
disp(['Re8 = ' num2str(Re8,'%10.10f')]);
disp(['De = ' num2str(De,'%10.10f')]);
disp(['Ca = ' num2str(Ca,'%10.10f')]);
disp(['LM = ' num2str(LAM,'%10.10f')]);
disp('--- Match ---');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% precomputation functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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