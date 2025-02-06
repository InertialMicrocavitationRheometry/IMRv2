close all
addpath('../../IMRv2/src/spectral/')
[t,R] = f_imrv2();%'r0',R0,'req',eqR);%,'mu',mu,'g',G);
hold on
plot(t,R)
