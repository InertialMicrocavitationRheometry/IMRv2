function varargout =  f_imrv2(varargin)
% IMR V2

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
% R - Bubble Radius 
% U - Bubble velocity 
% P - Internal bubble pressure
% T_Bubble - Temperature inside the bubble  
% T_Medium - Temperature outside the bubble  
% C - Vapor Concentration in the bubble
% Tm - Temperature in the medium 
% Dim - outputs variables in dimensional form
% Comp - 0 (ignores compressibility effects) or 1 (uses Keller- Miksis)
% Reduced - 0 utilizes full model or 1 uses Preston's reduced order model

%*************************************************************************
%   The various options can be found below.  Code may also run without any
%   inputs.
par = f_call_params(varargin{:});
params = cell2mat({par{1:end-1}});
tspan = cell2mat({par{end}});
i = 1;
% numerical settings 
radial          = params(i); i = i + 1;
bubtherm        = params(i); i = i + 1;
medtherm        = params(i); i = i + 1;
stress          = params(i); i = i + 1;
eps3            = params(i); i = i + 1;
masstrans       = params(i); i = i + 1;
% output options
dimensionalout  = params(i); i = i + 1;
progdisplay     = params(i); i = i + 1;
detail          = params(i); i = i + 1;
plotresult      = params(i); i = i + 1;
radiusonly      = params(i); i = i + 1;
vitalsreport    = params(i); i = i + 1;
displayonly     = params(i); i = i + 1;
technical       = params(i); i = i + 1;
% solver options
method          = params(i); i = i + 1;
spectral        = params(i); i = i + 1;
divisions       = params(i); i = i + 1;
% numerical parameters 
Nv              = params(i); i = i + 1;
Nt              = params(i); i = i + 1;
Mt              = params(i); i = i + 1;
Lv              = params(i); i = i + 1;
Lt              = params(i); i = i + 1;
% physical parameters%
% acoustic parameters
Cstar           = params(i); i = i + 1;
GAMa            = params(i); i = i + 1;
kappa           = params(i); i = i + 1;
nstate          = params(i); i = i + 1;
% dimensionless waveform parameters
tfin            = params(i); i = i + 1;
om              = params(i); i = i + 1;
ee              = params(i); i = i + 1;
tw              = params(i); i = i + 1;
dt              = params(i); i = i + 1;
mn              = params(i); i = i + 1;
wavetype        = params(i); i = i + 1;
pvarargin = {wavetype,om,ee,tw,dt,mn};
% dimensionless viscoelastic
We              = params(i); i = i + 1;
Re8             = params(i); i = i + 1;
DRe             = params(i); i = i + 1;
v_a             = params(i); i = i + 1;
v_nc            = params(i); i = i + 1;
Ca              = params(i); i = i + 1;
LAM             = params(i); i = i + 1;
De              = params(i); i = i + 1;
JdotA           = params(i); i = i + 1;
v_lambda_star   = params(i); i = i + 1;
iWe             = 1/We;
if Ca==-1; Ca=Inf; end
% dimensionless thermal 
Foh             = params(i); i = i + 1;
Br              = params(i); i = i + 1;
alpha           = params(i); i = i + 1;
beta            = params(i); i = i + 1;
chi             = params(i); i = i + 1;
iota            = params(i); i = i + 1;
% dimensionaless mass transfer 
Fom             = params(i); i = i + 1;
C0              = params(i); i = i + 1;
Rv_star         = params(i); i = i + 1;
Ra_star         = params(i); i = i + 1;
L_heat_star     = params(i); i = i + 1;
mv0             = params(i); i = i + 1;
ma0             = params(i); i = i + 1;
% dimensionless initial conditions
Rzero           = params(i); i = i + 1;
Uzero           = params(i); i = i + 1;
p0star          = params(i); i = i + 1;
P8              = params(i); i = i + 1;
T8              = params(i); i = i + 1;
Pv_star         = params(i); i = i + 1;
Req             = params(i);

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
sam = p0star - iWe + GAMa;
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
if stress == 1
    zeNO = 1; 
else 
    zeNO = 0; 
end
if spectral == 0, Nv = 1; end
if bubtherm == 0, Nt = -1; Mt = -1; qdot = []; end
if medtherm == 0, Mt = -1; end
ia = 4:(4+Nt);
ib = (5+Nt):(5+Nt+Mt);
ic = (6+Nt+Mt):(5+Nt+Mt+Nv);
id = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);

% initial condition assembly
init = [Rzero; Uzero; p0star; % radius, velocity, pressure
    zeros(Nt+1,1); % auxiliary temperature spectrum
    ones(Mt ~= -1); zeros(Mt,1); % medium temperature spectrum
    zeros(2*(Nv - 1)*(spectral == 1) + 2,1); % stress spectrum
    0]; % initial stress integral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SOLVER CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tspan = linspace(0,tfin,detail);
stepcount = 0;

if method == 15
    if divisions == 0
        options = odeset();
    else
        options = odeset('MaxStep',tfin/divisions,'RelTol',1e-6);
    end
    [t,X] = ode15s(@SVBDODE,tspan,init,options);
elseif method == 23
    if divisions == 0
        options = odeset();
    else
        options = odeset('MaxStep',tfin/divisions,'RelTol',1e-6);
    end
    [t,X] = ode23tb(@SVBDODE,tspan,init,options);
elseif method == 45
    if divisions == 0
        options = odeset('NonNegative',1,'AbsTol',1e-8,'RelTol',1e-8);
    else
        options = odeset('NonNegative',1,'MaxStep',tfin/divisions,'RelTol',1e-8);
    end
    [t,X] = ode45(@SVBDODE,tspan,init,options);
else
    if divisions == 0
        options = odeset('NonNegative',1);
    else
        options = odeset('NonNegative',1,'MaxStep',tfin/divisions,'RelTol',1e-6);
    end
    [t,X] = ode45(@SVBDODE,tspan,init,options);
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
    if bubtherm
        % extract auxiliary temperature
        SI = gA*X(ia);
        % auxiliary temperature derivatives
        dSI = gAPd*SI; % first order derivative
        ddSI = gAPdd*SI; % second order derivative
        % temperature and thermal diffusivity fields
        T = (alpha - 1 + sqrt(1+2*alpha*SI))/alpha;
        D = kapover*(alpha*T.^2 + (1-alpha)*T)/p;
        pVap = (f_pvsat(T(1)*T8)/P8);
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
        p = p0star*R^(-3*kappa);
        pdot = -3*kappa*U/R*p;
        pVap = Pv_star;        
    end

    J = 0; JdotX = 0; Z1dot = 0; Z2dot = 0;
    % stress equation
    if stress == 1 % Kelvin-Voigt with neo-Hookean elasticity
        % compute stress integral
        J = (4*(Req/R) + (Req/R)^4 - 5)/(2*Ca) - 4/Re8*U/R;
        JdotX = -2*U*(Req*(1/R)^2 + Req^4*(1/R)^5)/Ca + 4/Re8*U^2/R^2;
    elseif stress == 2 % linear Maxwell, linear Jeffreys, linear Zener
        % extract
        Z1 = X(ic);
        J = Z1/R^3 - 4*LAM/Re8*U/R;
        % stress integral derivative
        Z1dot = -Z1/De + 4*(LAM-1)/(Re8*De)*R^2*U - 4*(R^3-1)/(3*Ca*De);
        Z2dot = 0;
        JdotX = Z1dot/R^3 - 3*U/R^4*Z1 + 4*LAM/Re8*U^2/R^2;     
    elseif stress == 3 % upper-convected Maxwell, OldRoyd-B
        % extract stress sub-integrals
        Z1 = X(ic); Z2 = X(id);
        % compute new derivatives
        Z1dot = -(1/De - 2*U/R)*Z1 + 2*(LAM-1)/(Re8*De)*R^2*U;
        Z2dot = -(1/De + 1*U/R)*Z2 + 2*(LAM-1)/(Re8*De)*R^2*U;
        J = (Z1 + Z2)/R^3 - 4*LAM/Re8*U/R;
        JdotX = (Z1dot+Z2dot)/R^3 - 3*U/R^4*(Z1+Z2) + 4*LAM/Re8*U^2/R^2;
    elseif stress > 3 % Giesekus, PTT, or forced spectral         
        % extract stress spectrum
        c = X(ic); d = X(id);
        % inverse Chebyshev transforms and derivatives
        [trr,dtrr,t00,dt00] = stressdiff(c,d);
        % new spectral coefficient derivatives
        exptau = exp(ptt*Re8*De*(trr + 2*t00));            
        Z1dot = stresssolve(-(exptau/De + zeNO*4*U./(yV.^3*R) ...
            + eps3*Re8*trr).*trr ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dtrr ...
            - 4./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re8*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*4*LAM/Re8*(U/R)^2./yV.^6);
        Z2dot = stresssolve(-(exptau/De - zeNO*2*U./(yV.^3*R) ...
            + eps3*Re8*t00).*t00 ...
            + (1-ze).^2*U/(2*Lv*R).*(yV - zeNO./yV.^2).*dt00 ...
            + 2./yV.^3*((1-1/R^3)/(3*Ca) + U/(Re8*R) ...
            + LDR*(2*U^2/R^2 + Udot/R))/De - zeNO*10*LAM/Re8*(U/R)^2./yV.^6);
        % compute stress integral
        J = 2*sum(cdd.*(c-d));
        JdotX = 2*sum(cdd.*(Z1dot - Z2dot));           
    end
    
    % pressure waveform
    [pf8,pf8dot] = f_pinfinity(t,pvarargin{:});
  
    % bubble wall acceleration
    % Rayleigh-Plesset        
    if radial == 1
        Udot = (p + pVap - 1 - pf8 - iWe/R + J - 1.5*U^2)/R;
    % Keller-Miksis in enthalpy
    elseif radial == 2    
        hB = (sam/no)*(((p - iWe/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - iWe/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/Cstar)*(hB - pf8) - R/Cstar*pf8dot ...
            + R/Cstar*hH*(pdot + iWe*U./R.^2 + JdotX) ...
            - 1.5*(1 - U./(3*Cstar)).*U.^2)/((1 - U./Cstar).*R + JdotA*hH./Cstar);
    % Gilmore equation
	elseif radial == 3
        hB = sam/no*(((p - iWe/R + GAMa + J)/sam)^no - 1);
        hH = (sam/(p + pVap - iWe/R + GAMa + J))^(1/nstate);
        Udot = ((1 + U/Cstar)*(hB - pf8) - R/Cstar*pf8dot ...
            + R/Cstar*(hB + hH*(pdot + iWe*U/R^2 + JdotX)) ...
            - 1.5*(1 - U/(3*Cstar))*U^2) / ((1 - U/Cstar)*R + JdotA*hH/Cstar);
    % Keller-Miksis in pressure        
    elseif radial == 4  
        Udot = ((1+U./Cstar)*(p + pVap - 1 - pf8 - iWe./R + J) ...
            + R./Cstar.*(pdot + iWe.*U./R.^2 + JdotX - pf8dot) ...
            - 1.5.*(1-U./(3.*Cstar)).*U.^2)./((1-U./Cstar).*R + JdotA./Cstar);
    else 
        
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
a = X(:,ia)'; 
b = X(:,ib)'; 
c = X(:,ic)'; 
d = X(:,id)'; 
if masstrans == 1
    e = X(:,ie)';
end
I = X(:,end);
pA = zeros(size(t));
% Udot = R2dot(t,X);
for n = 1:length(t)
    pA(n) = f_pinfinity(t(n),pvarargin{:}); 
end 

% transform to real space
if spectral == 1
    trr = sCA*c; t00 = sCA*d;
else
    trr = c; t00 = d;
end
if bubtherm == 1
    T = (alpha-1+sqrt(1+2*alpha*gA*a))/alpha;
    if medtherm == 1
        TL = mA*b; 
    end
else
    T = R.^(-3*kappa);
end
if masstrans == 1
    C = gA*e;
end

if dimensionalout == 1
    
    % re-dimensionalize problem
    t = t*tc; R = R*Rref; U = U*uc; p = p*p0; pA = pA*p0; I = I*p0; T = T*T8;
    c = c*p0; d = d*p0; e = e*C0; C = C*C0; Udot = Udot*uc/tc;
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
        subplot(3,1,1);
        plot(t,R,'k','LineWidth',2); 
        hold('on'); 
        ylabel('$R$','Interpreter','Latex','FontSize',12);
        box on;
        axis([0 t(end) 0 (max(R) + min(R))]);
        set(gca,'TickLabelInterpreter','latex','FontSize',16)
        set(gcf,'color','w');
        subplot(3,1,3);
        hold on;
        box on;
        semilogy(t,abs(c(end,:)),'k-','LineWidth',2);
        semilogy(t,abs(d(end,:)),'b-','LineWidth',2); 
        xlabel('$t$','Interpreter','Latex','FontSize',12);
        ylabel('$c_P$, $d_P$','Interpreter','Latex','FontSize',12);
        set(gca, 'YScale', 'log');
        axis([0 t(end) 1e-20 1]);
        leg1 = legend('$c_P$','$d_P$','Location','NorthEast','FontSize',12);
        set(leg1,'Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex','FontSize',16)
        set(gcf,'color','w');
        if bubtherm == 1
            if medtherm == 0, b = zeros(size(a)); end
            subplot(3,1,2);
            hold on;
            box on;
            semilogy(t,abs(a(end,:)),'k-','LineWidth',2);
            semilogy(t,abs(b(end,:)),'b-','LineWidth',2); 
            ylabel('$a_N$, $b_M$, $e_N$','Interpreter','Latex','FontSize',12);
            set(gca, 'YScale', 'log');
            axis([0 t(end) 1e-20 1]);
            leg1 = legend('$a_N$','$b_M$','$e_N$','Location','NorthEast','FontSize',12);
            set(leg1,'Interpreter','latex');
        	set(gca,'TickLabelInterpreter','latex','FontSize',16)
            set(gcf,'color','w');
        end
    end
end

% assemble output
if displayonly == 1
    varargout = {};
else
    % standard outputs
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = U;
    varargout{4} = p;
    varargout{5} = trr;
    varargout{6} = t00;
    varargout{7} = I;
    varargout{8} = T;    
    if masstrans == 1
        varargout{9} = C;
    end
    if bubtherm == 1 && medtherm == 1
        varargout{10} = TL;
    else
        varargout{10} = ((T8 - 1)*dimensionalout + 1)*ones(divisions,1);
    end   
    % technical data
    if technical == 1
        varargout{11} = a; 
        varargout{12} = b;
        varargout{13} = c; 
        varargout{14} = d;
        varargout{15} = e;
        varargout{16} = [stepcount p0 alpha];
        varargout{17} = pA;
%         varargout{17} = Udot;
    end 
end

% convert run settings to strings
if radial == 1
    eqn = 'Rayleigh Plesset equation';
elseif enthalpy == 1
    eqn = 'Keller-Miksis in enthalpy';
else
    eqn = 'Keller-Miksis in pressure';
end
const = 'none';
if stress == 1
    if Ca == Inf
        const = 'Newtonian fluid';
    else
        const = 'neo-Hookean Kelvin-Voigt';
    end
elseif stress == 2
    if Ca ~= Inf && LAM == 0
        const = 'linear Zener';
    elseif Ca == Inf && LAM == 0
        const = 'linear Maxwell';
    elseif Ca == Inf && LAM ~= 0
        const = 'linear Jeffreys';
    else
        const = 'Kelvin-yangChurch series';
    end
elseif stress == 3
    if Ca ~= Inf && LAM == 0
        const = 'upper-convective Zener';
    elseif Ca == Inf && LAM == 0
        const = 'upper-convective Maxwell';
    elseif Ca == Inf && LAM ~= 0
        const = 'Oldroyd-B';
    end
elseif stress == 4 
    const = 'Phan-Thien-Tanner';
else 
    const = ['Giesekus(' num2str(eps3) ')'];
end

if bubtherm == 1
    if medtherm == 1, therm = 'full';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% functions called by solver %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
