clc;
clear;
close;

addpath('src');

% equation options
[tf,Rf,Uf] = m_imrv2_finitediff('Nt',100,'Mt',100);
[ts,Rs,Us] = m_imrv2_spectral('Nt',12,'Mt',12);