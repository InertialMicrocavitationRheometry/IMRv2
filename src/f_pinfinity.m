% file f_pinfinity.m
% brief contains function f_pinfinity

% brief This function the time-dependent external pressure in the
% surrounding material that drives the bubble oscillations. The function
% features wave types: histotripsy, Gaussian, and impulse
function [p8,p8dot] = f_pinfinity(t,vararg)
    
    % extracting input variables
    om = vararg(1);
    ee = vararg(2);
    tw = vararg(3);
    dt = vararg(4);
    mn = vararg(5);
    wave_type =  vararg(6);
    
    if wave_type < 0
        % experimentally-derived incoming waveform
        p8 = ppval(wave_poly,t);
        p8dot = ppval(wave_dpoly,t);
    elseif wave_type == 0
        % impulse wave
        [p8, p8dot] = impulse(t);
    elseif wave_type == 1
        % Gaussian wave
        [p8, p8dot] = gaussian(t);
    elseif wave_type == 2
        % histotripsy wave
        [p8, p8dot] = histo(t);
    elseif wave_type == 3
        % heaviside impulse
        [p8,p8dot] = heaviside_impulse(t);
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
    
    % Gaussian waveform
    function [p,pdot] = gaussian(t)
        p = - ee * exp(-((t-dt)^2) / tw^2);
        pdot = ee * (2*(t-dt)/tw^2) * exp(-((t-dt)^2) / tw^2);
    end
    
    % heaviside impulse
    function [p,pdot] = heaviside_impulse(t)
        p = -ee*(1-heaviside(t-tw));
        pdot = 0;
    end
    
end
