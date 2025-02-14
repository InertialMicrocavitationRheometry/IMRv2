% file f_pvsat.m
% brief contains function f_pvsat

% brief This function computes the saturated water vapor pressure per 
% (A. Preston 2004, Thesis)
function [ Pv ] = f_pvsat( T )
    Pv = 1.17e11*exp(-5200./(T)); 
end