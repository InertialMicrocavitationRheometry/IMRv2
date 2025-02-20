clc;
clear;
close;

addpath('src');

% equation options
[tf,Rf,Uf] = m_imrv2_finitediff('Nt',200,'Mt',200);
[ts,Rs,Us] = m_imrv2_spectral('Nt',12,'Mt',12);

figure(1)
hold on;
plot(tf,Rf,'^')
plot(ts,Rs,'o')
