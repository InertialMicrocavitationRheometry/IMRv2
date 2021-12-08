function varargout =  spectralkm2(varargin)
%GAMESETMATCH Solver for viscoelastic bubble dynamics
%   [t,R,p,pA,trr,t00,I,T,TL] = GeneralVisc;
%   [t,R,p,pA,trr,t00,I,T,TL] = GeneralVisc('prop1',x,'prop2',y,...);
%   [t,R,p,pA,trr,t00,I,T,TL,a,b,c,d,dat] = GeneralVisc('tech',1);
%   The various options can be found below.  Code may also run without any
%   inputs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default options and numerical parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display options
dimensionalout = 0; % output result in dimensional variables
progdisplay = 1; % display progress while code running
detail = 2000; % number of points in time to store result
plotresult = 1; % generate figure containing results
radiusonly = 0; % only produce R(t) curve
vitalsreport = 1; % display accuracy data

% output options
displayonly = 0; % do not generate output
technical = 0; % output technical data

% model for radial bubble dynamics, default is KME in pressure
rayleighplesset = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enthalpy = 0;

% thermal assumptions, default is no thermal assumptions
polytropic = 1;
cold = 0;
vapor = 0; % ignore vapor pressure

% constitutive model, default is UCM with linear elasticity
neoHook = 0;
voigt = 1;
linelas = 0;
liner = 0;
oldb = 0;
ptt = 0;
gies = 0;

% solver options
method = '45'; % ode45
spectral = 0; % force spectral collocation solution
divisions = 0; % minimum number of timesteps

% numerical parameters
Nv = 20;
Nt = 10;
Mt = 10;
Lv = 3;
Lt = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default physical parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% waveform parameters
TFin = 0.6e-6; % final time (s)
pA = 2e6; % pressure amplitude (Pa)
omega = 4e6*2*pi; % frequency (rad/s)
TW = 0; % gaussian width (s)
DT = 0; % delay (s)

% physical constants
GAM = 3049.13*1e5; % state equation parameter (Pa)
nstate = 7.15; % state equation parameter
p8 = 101325; % far-field pressure (Pa)
rho8 = 1060; % far-field density (kg/m^3)
c8 = sqrt(nstate*(p8 + GAM)/rho8); % far-field sound speed (m/s)
kappa = 1.4; % ratio of specific heats
S = 0.072; % surface tension (N/m)
T8 = 293.15; % far-field temperature (K)
cp = 4181; % specific heat of medium (J/(kg K))
pV = 2300;%2339.262; % vapor pressure of water (Pa)

% thermal conductivity and diffusivity
AT = 5.28e-5; % gas thermal conductivity slope parameter (W/(m K^2))
BT = 1.165e-2; % gas thermal conductivity shift parameter (W/(m K))
KL = 0.55; % medium thermal conductivity (W/(m K))
DL = 1.41e-7; % medium thermal diffusivity

% viscoelastic parameters
mu = 0.03; % viscosity (Pa s)
G = 20e3; % shear modulus (Pa) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda1 = 0.5e-6; % relaxation time (s)
lambda2 = lambda1; % retardation time (s)

% initial conditions
R0 = 1e-6; % initial radius (m)
U0 = 0; % initial velocity (m/s)
p0 = p8 + 2*S/R0 - pV*vapor; % initial pressure (Pa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overwrite defaults with options and dimensional inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that all inputs are matched
if mod(nargin,2) == 1, disp('Error: unmatched inputs'); return; end

% load inputs
misscount = 0; p0set = 0; c8set = 0;
if isempty(varargin) == 1
    disp('Using default settings');
else
    for n = 1:2:nargin
        if strcmpi(varargin{n},'p0') == 1, p0set = 1; end
        if strcmpi(varargin{n},'c8') == 1, c8set = 1; end
        switch lower(varargin{n})
            
            % thermal options
            case 'poly', polytropic = varargin{n+1} ~= 0;
            case 'cold', cold = varargin{n+1} ~= 0;
            case 'vapor', vapor = varargin{n+1} ~= 0;
            % equation for radial dynamics
            case 'rpe', rayleighplesset = varargin{n+1} ~= 0;
            case 'enth', enthalpy = varargin{n+1} ~= 0;
            case 'gil', gil = varargin{n+1} ~= 0;
                
            % constitutive model
            case 'neohook', neoHook = varargin{n+1} ~= 0;
            case 'voigt', voigt = varargin{n+1} ~= 0;
            case 'linelas', linelas = varargin{n+1} ~= 0;
            case 'liner', liner = varargin{n+1} ~= 0;
            case 'oldb', oldb = varargin{n+1} ~= 0;
            case 'ptt', ptt = varargin{n+1} ~= 0;
            case 'gies', gies = varargin{n+1};
                
            % display options
            case 'dimout', dimensionalout = varargin{n+1} ~= 0;
            case 'pdisp', progdisplay = varargin{n+1} ~= 0;
            case 'detail', detail = varargin{n+1};
            case 'plot', plotresult = varargin{n+1} ~= 0;
            case 'ronly', radiusonly = varargin{n+1} ~= 0;
            case 'vitals', vitalsreport = varargin{n+1} ~= 0;

            % output options
            case 'donly', displayonly = varargin{n+1} ~= 0;
            case 'tech', technical = varargin{n+1} ~= 0;

            % solver options
            case 'method', method = varargin{n+1};
            case 'spectral', spectral = varargin{n+1} ~= 0;
            case 'divisions', divisions = varargin{n+1};
                
            % numerical parameters
            case 'nv', Nv = varargin{n+1};
            case 'nt', Nt = varargin{n+1};
            case 'mt', Mt = varargin{n+1};
            case 'lv', Lv = varargin{n+1};
            case 'lt', Lt = varargin{n+1};
                
            % waveform parameters
            case 'tfin', TFin = varargin{n+1};
            case 'pa', pA = varargin{n+1};
            case 'omega', omega = varargin{n+1};
            case 'tw', TW = varargin{n+1};
            case 'dt', DT = varargin{n+1};
            case 'mn', mn = varargin{n+1};
                
            % physical constants
            case 'gam', GAM = varargin{n+1};
            case 'nstate', nstate = varargin{n+1};
            case 'p8', p8 = varargin{n+1};
            case 'rho8', rho8 = varargin{n+1};
            case 'c8', c8 = varargin{n+1};
            case 'kappa', kappa = varargin{n+1};
            case 'surf', S = varargin{n+1};
            case 't8', T8 = varargin{n+1};
            case 'pv', pV = varargin{n+1};
                
            % thermal conductivity and diffusivity
            case 'at', AT = varargin{n+1};
            case 'bt', BT = varargin{n+1};
            case 'kl', KL = varargin{n+1};
            case 'dl', DL = varargin{n+1};
                
            % viscoelastic parameters
            case 'mu', mu = varargin{n+1};
            case 'g', G = varargin{n+1};
            case 'lambda1', lambda1 = varargin{n+1};
            case 'lambda2', lambda2 = varargin{n+1};
                
            % initial conditions
            case 'r0', R0 = varargin{n+1};
            case 'rref', Rref = varargin{n+1};
            case 'u0', U0 = varargin{n+1};
            case 'p0', p0 = varargin{n+1};
            case 'a0', a0 = varargin{n+1};
            case 'b0', b0 = varargin{n+1};
                
            otherwise, misscount = misscount + 1;
        end
    end
    if p0set == 0, p0 = p8 + 2*S/R0 - pV*vapor; end
    if c8set == 0, c8 = sqrt(nstate*(p8 + GAM)/rho8); end
end

% check for physical viscoelastic parameters
if (lambda1 > mu/G && (voigt == 0 && neoHook == 0 && linelas == 0)) || abs(gies - 0.25) > 0.25
    disp('Error: nonphysical viscoelastic parameters');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nondimensionalize problem %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derived constants
uc = sqrt(p0/rho8); % characteristic velocity (m/s)
tc = Rref/uc; % characteristic time (s)
K8 = AT*T8 + BT; % far-field thermal conductivity (W/(m K))

% dimensionless state variables
C = c8/uc;
GAMa = GAM/p0;

% dimensionless waveform parameters
tfin = TFin/tc;
om = omega*tc;
ee = pA/p0;
tw = TW*tc;
dt = DT/tc;

% dimensionless vapor pressure
pVap = pV/p0*vapor;

% dimensionless numbers
We = 2*S/(p0*Rref); % Weber number
Re = p0*Rref/(mu*uc); % Reynolds number
Ca = p0/G; % Cauchy number
LAM = lambda2/lambda1;
De = lambda1*uc/Rref; % Deborah number
Fo = DL/(uc*Rref); % Fourier number
Br = uc^2/(T8*cp); % Brinkman number

% dimensionless thermal quantities
alpha = AT*T8/K8;
chi = T8*K8/(p0*Rref*uc);
iota = KL/(K8*Lt);

% dimensionless initial conditions
Rzero = R0/Rref;
Uzero = U0/uc;
pzero = p0*(R0/Rref)*(3*kappa)/p0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overwrite defaults with nondimensional inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) == 0
    for n = 1:2:nargin
        switch lower(varargin{n})
            
            % dimensionless state variables
            case 'cx', C = varargin{n+1};
            case 'gama', GAMa = varargin{n+1};
                
            % dimensionless waveform parameters
            case 'tfinx', tfin = varargin{n+1};
            case 'om', om = varargin{n+1};
            case 'ee', ee = varargin{n+1};
            case 'twx', tw = varargin{n+1};
            case 'dtx', dt = varargin{n+1};
                
            % dimensionless numbers
            case 'we', We = varargin{n+1};
            case 're', Re = varargin{n+1};
            case 'ca', Ca = varargin{n+1};
            case 'lam', LAM = varargin{n+1};
            case 'de', De = varargin{n+1};
            case 'fo', Fo = varargin{n+1};
                
            % dimensionless thermal quantities
            case 'alpha', alpha = varargin{n+1};
            case 'chi', chi = varargin{n+1};
            case 'iota', iota = varargin{n+1};
                
            % dimensionless initial conditions
            case 'rzero', Rzero = varargin{n+1};
            case 'uzero', Uzero = varargin{n+1};
            case 'pzero', pzero = varargin{n+1};
                
            otherwise, misscount = misscount + 1;
        end
    end
end

% check that all inputs were accounted for
if misscount ~= nargin/2
    disp(['Error: ' num2str(misscount-nargin/2) ' unrecognized input(s)']);
    return;
end

% check dimensionless viscoelastic parameters
if De > Ca/Re && (voigt == 0 && neoHook == 0 && linelas == 0)
    disp('Error: nonphysical viscoelastic parameters');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% final setting adjustments
if rayleighplesset == 1, enthalpy = 0; end
if enthalpy == 1, rayleighplesset = 0; end
if polytropic == 1, cold = 0; end
if cold == 1, polytropic = 0; end

if neoHook == 1
    [voigt,linelas,liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif voigt == 1
    [linelas,liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif linelas == 1
    [liner,ptt,gies,LAM] = deal(0);
    spectral = 0;
elseif liner == 1
    [voigt,ptt,gies] = deal(0);
elseif oldb == 1
    [voigt,liner,ptt,gies] = deal(0);
    Ca = Inf;
elseif ptt == 1
    [voigt,liner,gies,LAM] = deal(0);
    Ca = Inf; spectral = 1;
elseif gies ~= 0
    [voigt,liner,ptt,LAM] = deal(0);
    Ca = Inf; spectral = 1;
end

if voigt == 1 || neoHook == 1 || linelas == 1, JdotA = 4/Re;
elseif liner == 1 || oldb == 1, JdotA = 4*LAM/Re;
else JdotA = 0;
end
if spectral == 1, JdotA = 0; end

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
[sCA,sCI,sCAd,~,~,~] = dcdmtx(Nv);
sCA = sCA(2:end,2:end) - 1;
sCI = sCI(2:end,2:end);
sCAd = sCAd(2:end,2:end);

% precomputations
LDR = LAM*De/Re;
sam = 1 - We + GAMa;
no = (nstate-1)/nstate;
kapover = (kappa-1)/kappa;
yT = 2*Lt./(1+xi) - Lt + 1;
yV = 2*Lv./(1-ze) - Lv + 1;
nn = ((-1).^(0:Nt).*(0:Nt).^2)';
Udot = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if liner == 0, zeNO = 1; else zeNO = 0; end
if spectral == 0, Nv = 1; end
if polytropic == 1, Nt = -1; Mt = -1; qdot = []; end
if cold == 1, Mt = -1; end
ia = 4:(4+Nt);
ib = (5+Nt):(5+Nt+Mt);
ic = (6+Nt+Mt):(5+Nt+Mt+Nv);
id = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);

% initial condition assembly
init = [Rzero; Uzero; pzero; % radius, velocity, pressure
    zeros(Nt+1,1); % auxiliary temperature spectrum
    ones(Mt ~= -1); zeros(Mt,1); % medium temperature spectrum
    zeros(2*(Nv - 1)*(spectral == 1) + 2,1); % stress spectrum
    0]; % initial stress integral

% solver
tspan = linspace(0,tfin,detail);
stepcount = 0;

if strcmp(method,'15s')
    
    if divisions == 0, options = odeset();
    else options = odeset('MaxStep',tfin/divisions,'RelTol',1e-4);
    end
    
    [t,X] = ode15s(@SVBDODE,tspan,init,options);
    
elseif strcmp(method,'23s')
    
    if divisions == 0, options = odeset();
    else options = odeset('MaxStep',tfin/divisions);
    end
    
    [t,X] = ode23s(@SVBDODE,tspan,init,options);
    
elseif strcmp(method,'45')
    
    if divisions == 0, options = odeset('NonNegative',1);
    else options = odeset('NonNegative',1,'MaxStep',tfin/divisions,'RelTol',1e-12);
    end
    
    [t,X] = ode45(@SVBDODE,tspan,init,options);
    
else
    
    if divisions == 0, options = odeset('NonNegative',1);
    else options = odeset('NonNegative',1,'MaxStep',tfin/divisions);
    end
    
    [t,X] = ode45(@SVBDODE,tspan,init,options);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions called by solver %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if Nv < 650, s = sCI*x;
    else s = fctShift(x);
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
        p = R^(-3*kappa);
        pdot = -3*kappa*U/R*p;
        
    else
        
        % extract auxiliary temperature
        SI = gA*X(ia);
        
        % auxiliary temperature derivatives
        dSI = gAPd*SI; ddSI = gAPdd*SI;
        
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
            TLdot = (1+xi).^2/(Lt*R).*(Fo/R*((1+xi)/(2*Lt) - 1./yT) ...
                + U/2*(1./yT.^2 - yT)).*(mAPd*TL) ...
                + Fo/4*(1+xi).^4/(Lt^2*R^2).*(mAPdd*TL);
            
            % include viscous heating
%             if spectral == 1
%                 
%                 TLdot = TLdot - 2*Br*U./(R*yT.^3).*(ZZT*(X(ic) - X(id)));
%                 
%             elseif voigt == 1 || neoHook == 1 || linelas == 1
%                 
%                 TLdot = TLdot + 4*Br./yT.^6*(U/R*(1-1/R^3)/Ca + 3/Re*(U/R)^2);
%                 
%             end
            
            % enforce boundary condition and solve
            TLdot(end) = 0;
            qdot = [ones(1,Nt+1) -(alpha*(T(1)-1)+1)*ones(1,Mt+1); Q]...
                \[0; SIdot(2:end); TLdot(2:end); 0];
            
        end
    end
    
    % stress
    if neoHook == 1 % Kelvin-Voigt with neo-Hookean elasticity
        
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        
        % compute stress integral
        J = (4/R + 1/R^4 - 5)/(2*Ca) - 4/Re*U/R;
        JdotX = -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re*U^2/R^2;
        
    elseif voigt == 1 % Kelvin-Voigt (Yang-Church)
        
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        
        % compute stress integral
        J = -4/(3*Ca)*(1 - 1/R^3) - 4/Re*U/R;
        JdotX = -4/Ca*U/R^4 + 4/Re*U^2/R^2;
        
    elseif linelas == 1 % Linear elastic
        
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        
        % compute stress integral
        J = -2/Ca*(1 - 1/R^2) - 4/Re*U/R;
        JdotX = -4/Ca*U/R^3 + 4/Re*U^2/R^2;
        
    elseif spectral == 1 % Giesekus, PTT, or forced spectral
            
        % extract stress spectrum
        c = X(ic); d = X(id);

        % inverse Chebyshev transforms and derivatives
        [trr,dtrr,t00,dt00] = stressdiff(c,d);

        % new spectral coefficient derivatives
        exptau = exp(ptt*Re*De*(trr + 2*t00));            
        Z1dot = stresssolve(-(exptau/De + zeNO*4*U./(yV.^3*R) ...
            + gies*Re*trr).*trr ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dtrr ...
            - 4./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*4*LAM/Re*(U/R)^2./yV.^6);
        Z2dot = stresssolve(-(exptau/De - zeNO*2*U./(yV.^3*R) ...
            + gies*Re*t00).*t00 ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dt00 ...
            + 2./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*10*LAM/Re*(U/R)^2./yV.^6);

        % compute stress integral
        J = 2*sum(cdd.*(c-d));
        JdotX = 2*sum(cdd.*(Z1dot - Z2dot));
        
    elseif liner == 1 % linear Maxwell, linear Jeffreys, linear Zener
        
        % extract
        Z1 = X(ic);
        J = Z1/R^3 - 4*LAM/Re*U/R;

        % stress integral derivative
        Z1dot = -Z1/De + 4*(LAM-1)/(Re*De)*R^2*U - 4*(R^3-1)/(3*Ca*De);
        Z2dot = 0;
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re*U^2/R^2;     
        
        
    elseif oldb == 1 % upper-convected Maxwell, OldRoyd-B

        % extract stress sub-integrals
        Z1 = X(ic); Z2 = X(id);

        % compute new derivatives
        Z1dot = -(1/De - 2*U/R)*Z1 + 2*(LAM-1)/(Re*De)*R^2*U;
        Z2dot = -(1/De + 1*U/R)*Z2 + 2*(LAM-1)/(Re*De)*R^2*U;

        J = (Z1 + Z2)/R^3 - 4*LAM/Re*U/R;
        JdotX = (Z1dot+Z2dot)/R^3 - 3*U/R^4*(Z1+Z2) + 4*LAM/Re*U^2/R^2;
        
    end
    
    % bubble wall acceleration
    if rayleighplesset == 1
        
        % Rayleigh-Plesset
        Udot = (p + pVap - 1 + We - pf(t) - We/R + J - 1.5*U^2)/R;
       
        
    elseif enthalpy == 1
        
        % Keller-Miksis in enthalpy
        hB = sam/no*(((p - We/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - We/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/C)*(hB - pf(t)) - R/C*pfdot(t) ...
            + R/C*hH*(pdot + We*U/R^2 + JdotX) ...
            - 1.5*(1 - U/(3*C))*U^2)/((1 - U/C)*R + JdotA*hH/C);
        
	elseif gil == 1
        
        % TODO ADD GILMORE HERE
        hB = sam/no*(((p - We/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - We/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/C)*(hB - pf(t)) - R/C*pfdot(t) ...
            + R/C*hH*(pdot + We*U/R^2 + JdotX) ...
            - 1.5*(1 - U/(3*C))*U^2)/((1 - U/C)*R + JdotA*hH/C);
        
    else
        
        % Keller-Miksis in pressure
%         Udot = ((1+U/C)*(p + pVap - We/R + J - 1 + We - pf(t)) ...
%             + R/C*(pdot + We*U/R^2 + JdotX - pfdot(t)) ...
%             - 1.5*(1-U/(3*C))*U^2)/((1-U/C)*R + JdotA/C);
        % Keller-Miksis in pressure
        Udot = ((1+U./C)*(p + pVap - We./R + J - pf(t)) ...
            + R./C.*(pdot + We.*U./R.^2 + JdotX - pfdot(t)) ...
            - 1.5.*(1-U./(3.*C)).*U.^2)./((1-U./C).*R + JdotA./C);
        
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
Udot = R2dot(t,X);
for n = 1:length(t), pA(n) = pf(t(n)); end 

% transform to real space
if spectral == 1
    trr = sCA*c; t00 = sCA*d;
else
    trr = c; t00 = d;
end
if polytropic == 0
    T = (alpha-1+sqrt(1+2*alpha*gA*a))/alpha;
    if cold == 0, TL = mA*b; end
else
    T = R.^(-3*kappa);
end

if dimensionalout == 1
    
    % re-dimensionalize problem
    t = t*tc; R = R*Rref; U = U*uc; p = p*p0; pA = pA*p0; I = I*p0; T = T*T8;
    c = c*p0; d = d*p0; Udot = Udot*uc/tc;
    if spectral == 1, trr = trr*p0; t00 = t00*p0; end
    if polytropic == 0
        if cold == 0, TL = TL*T8; end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%
% output and display %
%%%%%%%%%%%%%%%%%%%%%%

% display figure
if plotresult == 1
    
    if vitalsreport == 0
        if radiusonly == 1
            plot(t,R); hold('on'); grid('on');
            axis([0 t(end) 0 (max(R) + min(R))]);
        else
            subplot(2 + spectral,1,1);
            plot(t,R); grid('on'); ylabel('R');
            axis([0 t(end) 0 (max(R) + min(R))]);
            subplot(2 + spectral,1,2);
            plot(t,I); ylabel('J');
            if spectral == 1
                subplot(3,1,3);
                plot(t,trr(end,:),t,t00(end,:));
                ylabel('\tau_{rr}|_R, \tau_{\theta\theta}|_R');
            end
            xlabel('t');
        end
    else
        
        subplot(3 - polytropic,1,1);
        plot(t,R); hold('on'); grid('on'); ylabel('R');
        axis([0 t(end) 0 (max(R) + min(R))]);
        subplot(3 - polytropic,1,2);
        semilogy(t,abs(c(end,:)),t,abs(d(end,:))); ylabel('c_P, d_P');
        axis([0 t(end) 1e-20 1]);
        if polytropic == 0
            if cold == 1, b = zeros(size(a)); end
            subplot(3,1,3);
            semilogy(t,abs(a(end,:)),t,abs(b(end,:))); ylabel('a_N, b_M');
            axis([0 t(end) 1e-20 1]);
        end
        xlabel('t');
    end
    
end

% assemble output
% if displayonly == 1, varargout = {};
% else
    
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
%     if technical == 1
        varargout{10} = a; varargout{11} = b;
        varargout{12} = c; varargout{13} = d;
        varargout{14} = [stepcount p0 alpha];
        varargout{15} = pA;
        varargout{16} = Udot;
%     end
    
% end

% convert run settings to strings
if rayleighplesset == 1, eqn = 'Rayleigh Plesset equation';
elseif enthalpy == 1, eqn = 'Keller-Miksis in enthalpy';
else eqn = 'Keller-Miksis in pressure';
end

if voigt == 0
    if neoHook == 1
        if Ca == Inf, const = 'Newtonian fluid';
        else const = 'neo-Hookean Voigt';
        end
    elseif linelas == 1
        if Ca == Inf, const = 'Newtonian fluid';
        else const = 'linear elastic Voigt';
        end
    elseif liner == 1
        if Ca ~= Inf && LAM == 0, const = 'linear Zener';
        elseif Ca == Inf && LAM == 0, const = 'linear Maxwell';
        elseif Ca == Inf && LAM ~= 0, const = 'linear Jeffreys';
        else const = 'Kelvin-Voigt series';
        end
    elseif liner == 0
        if ptt == 0 && gies == 0
            if Ca ~= Inf && LAM == 0, const = 'upper-convective Zener';
            elseif Ca == Inf && LAM == 0, const = 'upper-convective Maxwell';
            elseif Ca == Inf && LAM ~= 0, const = 'Oldroyd-B';
            end
        elseif Ca == Inf && LAM == 0
            if ptt == 1 && gies == 0, const = 'Phan-Thien-Tanner';
            elseif ptt == 0 && gies ~= 0
                const = ['Giesekus(' num2str(gies) ')'];
            end
        end
    end
elseif liner == 0 && gies == 0 && ptt == 0 && LAM == 0, const = 'Yang-Church Voigt';
end

if polytropic == 0
    if cold == 0, therm = 'full';
    else therm = 'cold-medium approximation';
    end
else therm = 'polytropic approximation';
end

if spectral == 1, solut = 'spectral method';
else solut = 'ODE formulation';
end

% display run settings
disp('--- Game settings ---');
disp(['Radial dynamics: ' eqn]);
disp(['Medium rheology: ' const]);
disp(['Thermal effects: ' therm]);
disp(['Solution method: ' solut]);
disp('--- Dimensionless numbers ---');
disp(['Re = ' num2str(Re,'%10.10f')]);
disp(['De = ' num2str(De,'%10.10f')]);
disp(['Ca = ' num2str(Ca,'%10.10f')]);
disp(['LM = ' num2str(LAM,'%10.10f')]);
disp('--- Match ---');

function dXdt = R2dot(t,X)

    if progdisplay == 1, disp(t/tfin); end
    
    % extract standard inputs
    R = X(:,1); U = X(:,2); p = X(:,3);

    % non-condensible gas pressure and temperature
    if polytropic == 1
        
        % polytropic approximation
        p = R.^(-3*kappa);
        pdot = -3*kappa*U./R.*p;
        
    else
        
        % extract auxiliary temperature
        SI = gA*X(ia);
        
        % auxiliary temperature derivatives
        dSI = gAPd*SI; ddSI = gAPdd*SI;
        
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
            TLdot = (1+xi).^2/(Lt*R).*(Fo/R*((1+xi)/(2*Lt) - 1./yT) ...
                + U/2*(1./yT.^2 - yT)).*(mAPd*TL) ...
                + Fo/4*(1+xi).^4/(Lt^2*R^2).*(mAPdd*TL);
            
            % enforce boundary condition and solve
            TLdot(end) = 0;
            qdot = [ones(1,Nt+1) -(alpha*(T(1)-1)+1)*ones(1,Mt+1); Q]...
                \[0; SIdot(2:end); TLdot(2:end); 0];
            
        end
    end
    
    % stress
    if neoHook == 1 % Kelvin-Voigt with neo-Hookean elasticity
        
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        
        % compute stress integral
        J = (4./R + 1./R.^4 - 5)./(2*Ca) - 4./Re.*U./R;
        JdotX = -2.*U.*(1./R.^2 + 1./R.^5)./Ca + 4./Re.*U.^2./R.^2;
        
    elseif voigt == 1 % Kelvin-Voigt (Yang-Church)
        
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        
        % compute stress integral
        J = -4/(3.*Ca)*(1 - 1./R.^3) - 4./Re.*U./R;
        JdotX = -4/Ca.*U./R.^4 + 4./Re*U.^2./R.^2;
        
    elseif linelas == 1 % Linear elastic
        
        % ignore stress sub-integrals
        Z1dot = 0; Z2dot = 0;
        
        % compute stress integral
        J = -2/Ca*(1 - 1/R.^2) - 4/Re*U./R;
        JdotX = -4/Ca*U./R.^3 + 4/Re*U.^2./R.^2;
        
    elseif spectral == 1 % Giesekus, PTT, or forced spectral
            
        % extract stress spectrum
        c = X(ic); d = X(id);

        % inverse Chebyshev transforms and derivatives
        [trr,dtrr,t00,dt00] = stressdiff(c,d);

        % new spectral coefficient derivatives
        exptau = exp(ptt*Re*De*(trr + 2*t00));            
        Z1dot = stresssolve(-(exptau/De + zeNO*4*U./(yV.^3*R) ...
            + gies*Re*trr).*trr ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dtrr ...
            - 4./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*4*LAM/Re*(U/R)^2./yV.^6);
        Z2dot = stresssolve(-(exptau/De - zeNO*2*U./(yV.^3*R) ...
            + gies*Re*t00).*t00 ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dt00 ...
            + 2./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*10*LAM/Re*(U/R)^2./yV.^6);

        % compute stress integral
        J = 2*sum(cdd.*(c-d));
        JdotX = 2*sum(cdd.*(Z1dot - Z2dot));
        
    elseif liner == 1 % linear Maxwell, linear Jeffreys, linear Zener
        
        % extract
        Z1 = X(ic);
        J = Z1/R^3 - 4*LAM/Re*U/R;

        % stress integral derivative
        Z1dot = -Z1/De + 4*(LAM-1)/(Re*De)*R^2*U - 4*(R^3-1)/(3*Ca*De);
        Z2dot = 0;
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re*U^2/R^2;     
        
        
    elseif oldb == 1 % upper-convected Maxwell, OldRoyd-B

        % extract stress sub-integrals
        Z1 = X(ic); Z2 = X(id);

        % compute new derivatives
        Z1dot = -(1/De - 2*U/R)*Z1 + 2*(LAM-1)/(Re*De)*R^2*U;
        Z2dot = -(1/De + 1*U/R)*Z2 + 2*(LAM-1)/(Re*De)*R^2*U;

        J = (Z1 + Z2)/R^3 - 4*LAM/Re*U/R;
        JdotX = (Z1dot+Z2dot)/R^3 - 3*U/R^4*(Z1+Z2) + 4*LAM/Re*U^2/R^2;
        
    end
    
    % bubble wall acceleration
    if rayleighplesset == 1
        
        % Rayleigh-Plesset
        Udot = (p + pVap - 1 + We - pf(t) - We/R + J - 1.5*U^2)/R;
       
        
    elseif enthalpy == 1
        
        % Keller-Miksis in enthalpy
        hB = sam./no.*(((p - We./R + GAMa + J)./sam).^no - 1);
        hH = (sam./(p + pVap - We./R + GAMa + J)).^(1/nstate);
        Udot = ((1 + U./C).*(hB - pf(t)) - R./C.*pfdot(t) ...
            + R./C.*hH.*(pdot + We*U./R.^2 + JdotX) ...
            - 1.5.*(1 - U./(3.*C)).*U.^2)./((1 - U./C).*R + JdotA.*hH./C);
        
	elseif gil == 1
        
        % Keller-Miksis in enthalpy
        hB = sam/no*(((p - We/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - We/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/C)*(hB - pf(t)) - R/C*pfdot(t) ...
            + R/C*hH*(pdot + We*U/R^2 + JdotX) ...
            - 1.5*(1 - U/(3*C))*U^2)/((1 - U/C)*R + JdotA*hH/C);
        
    else
        
        % Keller-Miksis in pressure
        Udot = ((1+U./C).*(p + pVap - We./R + J - pf(t)) ...
            + R./C.*(pdot + We.*U./R.^2 + JdotX - pfdot(t)) ...
            - 1.5.*(1-U./(3.*C)).*U.^2)./((1-U./C).*R + JdotA./C);
        
    end
    % output assembly
%     Jdot = JdotX - JdotA.*Udot./R;
    dXdt = Udot;
end




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