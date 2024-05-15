%% Clear everything
close all; clear all; clc;
%% Parameters
injections = 100;
mua = -3;
mub = 1;
mu = logspace(mua,mub,injections); 
Ga = -1;
Gb = 4;
G = P_inf./logspace(Ga,Gb,injections);  

%% IMRv2 all models 
mu = logspace(-3,-1,2);
G = logspace(-1,4,2);
model = {'neohook','voigt','linelas','liner','oldb'};
flag = [1,1,1,1,1];
addpath('../src/spectral')
IMR_data = struct();
for i = 1:length(model)
    for j = 1:length(mu)
        for k = 1:length(G)
            [t,R] = f_imrv2(model{i},flag(i),'mu',mu(j),'g',G(k));
            tvR_data = [t,R];
            IMR_data(j,k,i).tvR = tvR_data;
            mu(j)
            G(k)
        end
    end
end
save('IMR_data.mat','IMR_data')

%% Scratch paper


[t,R] = f_imrv2('voigt',1)






