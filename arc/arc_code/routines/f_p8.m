function [p8,dp8] = f_p8(t,f,deltap,force,po)
% FORCING PRESSURE CALCULATION
 if strcmp('sine',force)==1
         p = deltap*(sin(2*pi*f*(t)));
 elseif strcmp('mono',force)==1
     p = deltap;
 end     
 p8 = p + po;

 if strcmp('sine',force)==1
         dp8 = deltap*2*pi*f*cos(2*pi*f*(t));
 elseif strcmp('mono',force)==1
     dp8 = 0;
 end 

end

