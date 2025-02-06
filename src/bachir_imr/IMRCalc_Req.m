function [R_eq,P_eq,C_eq] = IMRCalc_Req(R0,Tgrad,Cgrad,Pa,G,G1,mu)

%***********************************************************************
% NOTE: 
% Used to estimate the equilibrium radii 
% It assumed that at equilibrium the medium is stress free, so
% R_eq is needed to set up the Elastic force on the bubble f(R_eq/R) 
%***********************************************************************

%***************************************
% Load Parameters : 
Pmt = IMRcall_parameters(R0,G,G1,mu); % Calls parameters script 
k = Pmt(1);  We = Pmt(7);  Rv_star = Pmt(11);  Ra_star = Pmt(12); 
P_inf = Pmt(19); T_inf = Pmt(20);

%***************************************

Pv = Pvsat(1*T_inf)/P_inf; Pa = Pa/P_inf; ma0 = Pa/Ra_star;
MTotal0 = Pa/Ra_star + Pv/Rv_star;


if Cgrad == 0 && Tgrad == 0
    
    fun = @(x) Pa *(1/x)^(3*k)-(1)+Pv-1/(We*x);

    exp2 = 1;
    x = fzero(fun,exp2);

    while (isnan(x))

        exp2 = exp2*1.01;
        x = fzero(fun,exp2);

    end

    R_eq = x;
    P_eq = Pa*(1/R_eq)^(3*k)+Pv; 
    theta = Rv_star/Ra_star*(P_eq/Pv-1);
    C_eq = 1/(1+theta);
    
elseif Cgrad == 0 && Tgrad == 1
    
    % Note: For very small initial mass content this case has a hard time
    % finding roots. One fix is to just set Cgrad to 1. That always
    % converges 
    %fun = @(x) Pa*(1/x)^(3)-(1)+Pv-1/(We*x);
    fun = @(x) Pa + (-1+Pv)*x^3 -x^2/(We);
    exp2=1;
    x = fzero(fun,exp2);
    while (isnan(x))
        exp2 = exp2*1.01;
        x = fzero(fun,exp2);
    end
    R_eq = x;
    P_eq = Pa*(1/ R_eq )^(3)+Pv;
    theta = Rv_star/Ra_star*(P_eq/Pv-1);
    C_eq = 1/(1+theta);
        
elseif Cgrad == 1 && Tgrad == 1
    
   % Full model:    
    
    fun = @(x) Pv*(1+(ma0/x)*(Ra_star/Rv_star))-1-...
    1/We*(Pv/(Rv_star*x))^(1/3) ; % parameterized function

    
    exp2=MTotal0;
    x = fzero(fun,exp2,optimset('display','off'));

    while (isnan(x))

        exp2 = exp2/1.11;
        x = fzero(fun,exp2,optimset('display','off'));

    end
    MVE = x;
    R_eq  =(Rv_star*MVE/Pv)^(1/3);
    P_eq = Pv*(1+ma0/MVE*Ra_star/Rv_star);
    theta = Rv_star/Ra_star*(P_eq/Pv-1);
    C_eq = 1/(1+theta); 



end 

end