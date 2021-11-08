function [mu8,Dmu] = f_nonNewtonian_Re(vmaterial)
%F_NONNEWTONIAN Outputs the Reynolds number that is dynamically changes
% with the shear rate. Note: Re = P_inf*R0/(m8*Uc). Units are in Pascal
% seconds.

    if strcmp('water',vmaterial)==1
        mu8 = 8.3283e-4;
        muo = 8.3283e-4;
    elseif strcmp('blood', vmaterial) == 1
        mu8 = 0.00345; 
        muo = 0.056; 
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
    end
    
    Dmu = muo-mu8;
    
end