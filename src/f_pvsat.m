function [ Pv ] = f_pvsat( T )
%Calculates the saturated vapor pressure using temperature 
% Per (A. Preston 2004, Thesis)

Pv = 1.17e11*exp(-5200./(T)); 

end