function [p8,p8dot] = f_pinfinity(t,vararg)
%F_PINFINITY Summary of this function goes here
%   Detailed explanation goes here

om =        vararg(1);
ee =        vararg(2);
tw =        vararg(3);
dt =        vararg(4);
mn =        vararg(5);
wavetype =  vararg(6);

switch wavetype
    case 0
        [p8, p8dot] = histo(t);
    case 1
        [p8, p8dot] = gaussian(t);
    case 2
        [p8, p8dot] = impulse(t);
end

% histotripsy waveform
function [p,pdot] = histo(t)
    if t < dt - pi/om
        p = 0;
    elseif t > dt + pi/om
        p = 0;
    else
        p = ee*(0.5 + 0.5*cos(om*(t - dt))).^mn;
    end
    if t < dt - pi/om
        pdot = 0;
    elseif t > dt + pi/om
        pdot = 0;
    else
        pdot = -ee*mn*(0.5+0.5*cos(om*(t-dt))).^(mn-1)*0.5*om.*sin(om*(t-dt));
    end    
end

% impulse waveform
function [p,pdot] = impulse(~)
    p = ee;
    pdot = 0;
end

end