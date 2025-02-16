function [mu8,Dmu,a,nc,lambda,vmat] = f_nonNewtonian_Re(vmaterial)
    %F_NONNEWTONIAN Outputs the Reynolds number that is dynamically changes
    % with the shear rate. Note: Re = P8*R0/(m8*Uc). Units are in Pascal
    % seconds.
    a = 0;
    nc = 0;
    lambda = 0;
    if strcmp('water',vmaterial)==1
        mu8 = 8.3283e-4;
        muo = 8.3283e-4;
        vmat = 1;
    elseif strcmp('blood_infinity', vmaterial) == 1
        mu8 = 0.00345;
        muo = mu8;
        vmat = 2;
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
        error('INPUT ERROR: No viscosity model specified in f_call_parameters, exiting');
    end
    Dmu = muo-mu8;
end
