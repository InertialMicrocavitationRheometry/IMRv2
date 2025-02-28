clc;
clear;
close;

addpath('src');

% equation options
[tf,Rf,Uf] = m_imr_finitediff('Nt',100,'Mt',100);
[ts,Rs,Us] = m_imr_spectral('Nt',12,'Mt',12);