% file f_pinfinity.m
% brief contains function f_pinfinity

% brief This function the time-dependent external pressure in the
% surrounding material that drives the bubble oscillations. The function
% features wave types: histotripsy, Gaussian, and impulse
function [p8,p8dot] = f_pinfinity(t,vararg)
    %F_PINFINITY Summary of this function goes here
    %   Detailed explanation goes here
    
    om = vararg(1);
    ee = vararg(2);
    tw = vararg(3);
    dt = vararg(4);
    mn = vararg(5);
    wave_type =  vararg(6);
    
    switch wave_type
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

%TODO
% if (Pext_type == 'HN')
%     in = load('./data/workspace.mat','pp_HN');
%     pp_HN = in.pp_HN;
%     in = load('./data/workspace.mat','dp_HN');
%     dp_HN = in.dp_HN;
% elseif (Pext_type == 'MN')
%     in = load('./data/workspace.mat','pp_MN');
%     pp_MN = in.pp_MN;
%     in = load('./data/workspace.mat','dp_MN');
%     dp_MN = in.dp_MN;
% elseif (Pext_type == 'ML')
%     in = load('./data/workspace.mat','pp_ML');
%     pp_ML = in.pp_ML;
%     in = load('./data/workspace.mat','dp_ML');
%     dp_ML = in.dp_ML;
% else
% end

%TODO
% % set external pressure
% if (Pext_type == 'sn')
%     Pext =  -Pext_Amp_Freq(1)/P_inf*sin(2*pi*Pext_Amp_Freq(2)*t*t0) ;
%     P_ext_prime = -2*pi*Pext_Amp_Freq(2)*t0*Pext_Amp_Freq(1)/P_inf...
    %         *cos(2*pi*Pext_Amp_Freq(2)*t*t0) ;
%
% elseif (Pext_type == 'RC')
%
%     Pext = Pext_Amp_Freq(1)/P_inf ;
%            P_ext_prime = 0;
%
% elseif (Pext_type == 'GS')
%     a = 10;
%     tw = 1;
%     Pext = -Pext_Amp_Freq(1)*(exp(-((t-a)/tw)^2))/P_inf;
%     P_ext_prime = 2*Pext_Amp_Freq(1)/P_inf.*(exp(-((t-a)./tw).^2)).*(t*t0-a)./tw^2;
%
% elseif (Pext_type == 'RG')
%
%     Pext = -Pext_Amp_Freq(1)/P_inf ;
%     P_ext_prime = 0;
%
% elseif (Pext_type == 'ip')
%
%     Pext = -Pext_Amp_Freq(1)/P_inf*...
    %     (1-heaviside(t-Pext_Amp_Freq(2)/t0)) ;
%     P_ext_prime = 0;
%
%  elseif (Pext_type == 'IC')
%
%     Pext = 0;
%     P_ext_prime = 0;
%
% elseif (Pext_type == 'HN')
%     Pext = ppval(pp_HN,t);
%     P_ext_prime = ppval(dp_HN,t);
%
% elseif (Pext_type == 'MN')
%     Pext = ppval(pp_MN,t);
%     P_ext_prime = ppval(dp_MN,t);
%
% elseif (Pext_type == 'ML')
%     Pext = ppval(pp_ML,t);
%     P_ext_prime = ppval(dp_ML,t);
%
% end
