clc;
clear;
close;

addpath('src');

% equation options
[tf,Rf,Uf] = m_imrv2_finitediff('Nt',20,'Mt',20);
[ts,Rs,Us] = m_imrv2_spectral('Nt',20,'Mt',20);