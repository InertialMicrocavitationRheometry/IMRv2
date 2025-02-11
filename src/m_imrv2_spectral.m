% file m_imrv2_spectral.m
% brief contains module m_imrv2_spectral

% brief This module features a Chebyshev spectral collocation solver of the
% PDEs involving thermal transport and viscoelasticity to solve
% Rayleigh-Plesset equations
function varargout =  m_imrv2_spectral(varargin)

% problem Initialization
[eqns_opts, solve_opts, init_opts, tspan_opts, out_opts, acos_opts,... 
    wave_opts, sigma_opts, thermal_opts, mass_opts]...
    = f_call_params(varargin{:});

% equations settings 
radial          = eqns_opts(1);  bubtherm        = eqns_opts(2); 
medtherm        = eqns_opts(3);  stress          = eqns_opts(4); 
eps3            = eqns_opts(5);  vapor           = eqns_opts(6);
masstrans       = eqns_opts(7);  
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
Req             = init_opts(7);  alphax          = init_opts(8);

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

% pre_process code

% collocation point construction
y = cos(pi*(0:Nt)'/(2*Nt));
xi = cos(pi*(0:Mt)'/Mt);
ze = cos(pi*(1:Nv)'/Nv);
% collocation matrix construction
[gA,gAI,~,~,gAPd,gAPdd] = dcdmtxe(Nt);
[mA,~,~,~,mAPd,mAPdd] = dcdmtx(Mt);
[gC,gCI,~,~,~,~] = dcdmtxe(Nt);
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
sam = 1 - Pv_star + GAMa; 
no = (nstate-1)/nstate;
kapover = (kappa-1)/kappa;
yT = 2*Lt./(1+xi) - Lt + 1;
yV = 2*Lv./(1-ze) - Lv + 1;
nn = ((-1).^(0:Nt).*(0:Nt).^2)';
nn = sparse(nn);
Udot = 0; 

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
if stress == 1
    zeNO = 0; 
else 
    zeNO = 1; 
end
if spectral == 0, Nv = 1; end
if bubtherm == 0, Nt = -1; Mt = -1; qdot = []; end
if medtherm == 0, Mt = -1; end
if masstrans == 0, Nm = -1; end
ia = 4:(4+Nt);
ib = (5+Nt):(5+Nt+Mt);
ic = (6+Nt+Mt):(5+Nt+Mt+Nv);
id = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);
ie = (6+Nt+Mt+2*Nv):(5+Nt+Mt+2*Nv+Nm);

% initial condition assembly

% radius, velocity, pressure, auxiliary temperature spectrum,
% medium temperature spectrum, stress spectrum,
% initial stress integral
Tau0 = zeros(Nt+1,1);
Tm0 = ones(Mt ~= -1);
Tm1 = zeros(Mt,1);
% TODO ADD THE MASS Transfer structure here

Sp = zeros(2*(Nv - 1)*(spectral == 1) + 2,1);  
init = [Rzero; Uzero; p0star; Tau0; Tm0; Tm1; Sp; 0]; 

% solver 
f_display(radial, bubtherm, medtherm, masstrans, stress, spectral, eps3, Re8, De, Ca, LAM);
stepcount = 0;
bubble = @SVBDODE;
[t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);

% post processing

% extract result
R = X(:,1); 
U = X(:,2); 
p = X(:,3); 
% extracting the Chebyshev coefficients
Z1 = X(:,ic); 
Z2 = X(:,id); 
a = X(:,ia)';
b = X(:,ib)'; 
c = X(:,ic)'; 
d = X(:,id)'; 
e = X(:,ie)';
I = X(:,ie+1);
if bubtherm
    T = (alpha-1+sqrt(1+2*alpha*gA*a))/alpha;
    if medtherm
        TL = mA*b; 
    end
else
    T = R.^(-3*kappa);
end
if masstrans
    C = gC*e;
end

pA = zeros(size(t));
for n = 1:length(t)
    pA(n) = f_pinfinity(t(n),pvarargin); 
end 

% transform to real space
if spectral == 1
    trr = sCA*c; t00 = sCA*d;
else
    trr = c; t00 = d;
end
% dimensionalization
if dimensionalout == 1
    % re-dimensionalize problem
    t = t*tc; 
    R = R*Rref; 
    U = U*uc; 
    p = p*p0; 
    T = T*T8;
    %pA = pA*p0; 
    I = I*p0; 
    c = c*p0; 
    d = d*p0; 
    e = e*C0; 
    C = C*C0; 
    Udot = Udot*uc/tc;
    if spectral == 1
        trr = trr*p0; 
        t00 = t00*p0; 
    end
    if bubtherm == 1
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
varargout{8} = trr;
varargout{9} = t00;
varargout{10} = I;
if stress == 2
    varargout{11} = Z1; % Z1
    varargout{12} = Z2; % Z2
    varargout{13} = (Z1)./R.^3 - 4*LAM/Re8.*U./R; % J
elseif stress == 3
    varargout{11} = Z1; % Z1
    varargout{12} = Z2; % Z2
    varargout{13} = (Z1 + Z2)./R.^3 - 4*LAM./Re8.*U./R; % J
end

% solver function
function dXdt = SVBDODE(t,X)
    stepcount = stepcount + 1;
    if progdisplay == 1, disp(t/tfin); end
    
    % extract standard inputs
    R = X(1); 
    U = X(2); 
    p = X(3); 
    qdot = []; 
    
    % non-condensible gas pressure and temperature
    if bubtherm
        % extract auxiliary temperature
        SI = gA*X(ia);
        % auxiliary temperature derivatives
        dSI = gAPd*SI; % first order derivative
        ddSI = gAPdd*SI; % second order derivative
        % temperature and thermal diffusivity fields
        T = (alpha - 1 + sqrt(1+2*alpha*SI))/alpha;
        D = kapover*(alpha*T.^2 + (1-alpha)*T)/p;
        pVap = vapor*(f_pvsat(T(1)*T8)/P8);

        % new pressure and auxiliary temperature derivative
        pdot = 3/R*((kappa-1)*chi/R*dSI(1) - kappa*p*U);
        
        SIdot = pdot*D + chi/R^2*(2*D./y - kapover/p*(dSI - dSI(1)*y)).*dSI ...
            + chi*D/R^2.*ddSI;
        SIdot(end) = pdot*D(end) - chi/R^2*(8*D(end)*sum(nn.*X(ia)) ...
            + kapover/p*dSI(end)^2) + chi*D(end)/R^2.*ddSI(end);
        if medtherm % warm-liquid
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
            elseif stress == 1
                TLdot = TLdot + 4*Br./yT.^6*(U/R*(1-1/R^3)/Ca + 3/Re8*(U/R)^2); 
            end
            % enforce boundary condition and solve
            TLdot(end) = 0;
            qdot = [ones(1,Nt+1) -(alpha*(T(1)-1)+1)*ones(1,Mt+1); Q]...
                \[0; SIdot(2:end); TLdot(2:end); 0];       
        else % cold-liquid approximation
            % solve auxiliary temperature with boundary condition
            qdot = gAI*[0; SIdot(2:end)];            
        end
    else 
        % polytropic approximation
        p = (p0star-Pv_star)*R^(-3*kappa);
        pdot = -3*kappa*U/R*p;
        pVap = Pv_star;        
    end
    
    % stress equation
    Z1dot = 0; Z2dot = 0;

        %TODO Need to add non-Newtonian behavior to JdotX 
        %((1-U/C_star)*R + ...
        %  4/Re8/C_star - 6*ddintfnu*iDRe/C_star);

    % no stress
    if stress == 0
        J = 0;
        JdotX = 0;

    % compute stress integral
    elseif stress == 1 
        % Kelvin-Voigt with neo-Hookean elasticity
        J = (4*(Req/R) + (Req/R)^4 - 5)/(2*Ca) - 4/Re8*U/R;
        JdotX = -2*U*(Req*(1/R)^2 + Req^4*(1/R)^5)/Ca + 4/Re8*U^2/R^2;

    elseif stress == 2
        % quadratic Kelvin-Voigt with neo-Hookean elasticity
        J = (3*alphax-1)*(5 - (Req/R)^4 - 4*(Req/R))/(2*Ca) - 4/Re8*U/R + ...
        (2*alphax/Ca)*(27/40 + (1/8)*(Req/R)^8 + (1/5)*(Req/R)^5 + (1/2)*(Req/R)^2 - ...
        2*R/Req);
        JdotX = ((3*alphax-1)/(2*Ca))*((4*Req^4*U/R^5) + (4*Req*U/R^2)) + ...
        4*(U^2)/(Re8*R^2) - (2*alphax/Ca)*(2*U/Req + Req^8*U/R^9 + ...
        Req^5*U/R^6 + Req^2*U/R^3);

    % Giesekus, PTT, or forced spectral             
    elseif spectral
        % extract stress spectrum
        c = X(ic); d = X(id); 
        % inverse Chebyshev transforms and derivatives
        [trr,dtrr,t00,dt00] = stressdiff(c,d);
        % new spectral coefficient derivatives
        exptau = exp(ptt*Re8*De*(trr + 2*t00));            
        Z1dot = stresssolve(-(exptau/De + zeNO*4*U./(yV.^3*R) ...
            + eps3*Re8*trr).*trr ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dtrr ...
            - 4./yV.^3*((1-(Req/R)^3)/(3*Ca) + U/(Re8*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*4*LAM/Re8*(U/R)^2./yV.^6);
        Z2dot = stresssolve(-(exptau/De - zeNO*2*U./(yV.^3*R) ...
            + eps3*Re8*t00).*t00 ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dt00 ...
            + 2./yV.^3*((1-(Req/R)^3)/(3*Ca) + U/(Re8*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*10*LAM/Re8*(U/R)^2./yV.^6);
        % compute stress integral
        J = 2*sum(cdd.*(c-d));
        JdotX = 2*sum(cdd.*(Z1dot - Z2dot));   

    % linear Maxwell, linear Jeffreys, linear Zener        
    elseif stress == 3 
        % extract
        Z1 = X(ic);
        J = Z1/R^3 - 4*LAM/Re8*U/R;
        Ze = R^3*((3*alphax-1)*(5 - (Req/R)^4 - 4*(Req/R))/(2*Ca) + ...
        (2*alphax/Ca)*(27/40 + (1/8)*(Req/R)^8 + (1/5)*(Req/R)^5 + (1/2)*(Req/R)^2 - ...
        2*R/Req));
        ZdotSqNH = 3*R^2*U*((3*alphax-1)*(5 - (Req/R)^4 - 4*(Req/R))/(2*Ca) + ...
        (2*alphax/Ca)*(27/40 + (1/8)*(Req/R)^8 + (1/5)*(Req/R)^5 + (1/2)*(Req/R)^2 - ...
        2*R/Req)) + R^3*(((3*alphax-1)/(2*Ca))*((4*Req^4*U/R^5) + (4*Req*U/R^2)) - ...
        (2*alphax/Ca)*(2*U/Req + Req^8*U/R^9 + ...
        Req^5*U/R^6 + Req^2*U/R^3)); % ddt(R^3 S_qKV)
        % ZdotNH = -1/(2*Ca)*(3*R^2*U*(5-(Req/R)^4-4*Req/R)+ ...
            % R^2*U*(4*(Req/R)^5+4*Req/R));
        % ZdotYC = - 4*(R^3-Req^3)/(3*Ca*De);
        % stress integral derivative
        Z1dot = -(Z1-Ze)/De + ZdotSqNH + (3*U/R)*(Z1-Ze) + 4*(LAM-1)/(Re8*De)*R^2*U ;
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re8*U^2/R^2;

    % upper-convected Maxwell, OldRoyd-B
    elseif stress == 4 
        % extract stress sub-integrals
        Z1 = X(ic); Z2 = X(id);
        % compute new derivatives
        Z1dot = -(1/De - 2*U/R)*Z1 + 2*(LAM-1)/(Re8*De)*R^2*U;
        Z2dot = -(1/De + 1*U/R)*Z2 + 2*(LAM-1)/(Re8*De)*R^2*U;
        J = (Z1 + Z2)/R^3 - 4*LAM/Re8*U/R;
        JdotX = (Z1dot+Z2dot)/R^3 - 3*U/R^4*(Z1+Z2) + 4*LAM/Re8*U^2/R^2;        
    else
        error('stress setting is not available');
    end

    % pressure waveform
    [pf8,pf8dot] = f_pinfinity(t,pvarargin);
    
    % bubble wall acceleration
    [Udot] = f_radial_eq(radial, p, pdot, pVap, pf8, pf8dot, iWe, R, U, ...
        J, JdotX, Cstar, sam, no, GAMa, nstate, JdotA );

    % stress integral rate
    Jdot = JdotX - JdotA*Udot/R;

    % output assembly
    dXdt = [U; Udot; pdot; qdot; Z1dot; Z2dot; Jdot];

end
disp('--- COMPLETED SIMULATION ---');  

% functions called by solver 

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

% function Cw= CW(Tw,P)
%   % Calculates the concentration at the bubble wall 
%   %Function of P and temp at the wall 
%   thetha = Rv_star/Ra_star*(P./(f_pvsat(Tw*T8)/P8) -1);
%   Cw = 1./(1+thetha); 
% end

end

% precomputation functions 
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
