function [vecoutput] = f_multirun()
close all; clear all; clc;
%consider neoHook and linear elastic 
%for Newtonian and elastic materials
%viscosity will range from zero to 10^1
%shear modulus will range from zero to 10^3

addpath('src/spectral');

%starting with linear elasticity
%linear elasticity = 1 ; neoHook = 0
muvec = [0,logspace(-6,1,10)];
Gvec = [0,logspace(-3,3,10)];

% for i=1:length(muvec)
%     for j = 1:length(Gvec)
%         mu = muvec(i);
%         G = Gvec(j);
%         foldername = '../bdata/linearelasticity';
%         mkdir(foldername);
%         filename = strcat(foldername,'/data_',num2str(mu),'_',num2str(G),'.mat');
%         [t,R,p,Rdot,trr,t00,~,T,C,TL] = f_imrv2('linelas',1,'neoHook',0,...
%          'mu',mu,'g',G);
%         save(filename);
%     end
% end
%neoHook = 1; linear elasticity =0
foldername = '../3dasm_data/sim_data/neoHookean';
mkdir(foldername);
for i=1:length(muvec)
    for j = 1:length(Gvec)
        mu = muvec(i);
        G = Gvec(j);
        filename = strcat(foldername,'/data_',num2str(mu),'_',num2str(G),'.csv');
        [t,R,U,P]=f_imrv2('linelas',0,'neoHook',1,'mu',mu,'g',G);
        
        x = 2*(t./(t(end)-t(1)))- (t(end)+t(1))/(t(end)-t(1));
        y = 2*(R./(max(R)-min(R)))-(max(R)+min(R))/(max(R)-min(R));
        a = [x,y];
        writematrix(a, filename);
        plot(x,y,'s')
    end
end
%% This is for 10% PA/0.06% BIS experimental data
clear all; close all; clc;
%addpath("../3dasm_data/exp_data")
load("PA_10%_0.06%_completed.mat");
%extract quantities of interest only for the first file
a = []; b=[];
for kk = 1:37
    file = expts(kk);
    t = file.t_norm;
    R = file.R_norm;
    x = 2*(t./(t(end)-t(1)))- (t(end)+t(1))/(t(end)-t(1));
    y = 2*(R./(max(R)-min(R)))-(max(R)+min(R))/(max(R)-min(R));
    a = [a;x;y];
    
    t0 = file.t0;
    R0 = file.R0;
    b = [b;t0;R0];
end
a = a';
t_all = a(:,1:2:end-1);
t_avg = mean(t_all,2);
R_all = a(:,2:2:end);
R_avg = mean(R_all,2);
c = [t_avg,R_avg];

d = d';
t0_all = b(:,1:2:end-1);
R0_all = b(:,2:2:end);
t0_avg = mean(t0_all,2);
R0_avg = mean(R0_all,2);

%filename = strcat('file_,num2str(kk),'.csv');
writematrix(c,"Experimental_Data_10PA_.06BIS.csv")
R_std = std(R_all,1,2);

%fill([t_avg,flip(t_avg)],[R_avg+R_std,flip(R_avg-R_std)],'g')
%patch(t_avg,[R_avg+R_std,flip(R_avg-R_std)],'g')
figure(1)
hold on
plot(t_avg,R_avg,'b')
plot(t_avg,R_avg-R_std,'k')
plot(t_avg,R_avg+R_std,'r')
%fill([t_avg,flip(t_avg)],[R_avg+R_std,flip(R_avg-R_std)],'g')


%% This is for one simulated data with the best fit G and mu
G_best = 1.39e+04;
mu_best = 0.1070;

[t,R,U,P]=f_imrv2('linelas',0,'neoHook',1,'mu',mu_best,'g',G_best,'t0',t0_avg,'R0',R0_avg);
x = 2*(t./(t(end)-t(1)))- (t(end)+t(1))/(t(end)-t(1));
y = 2*(R./(max(R)-min(R)))-(max(R)+min(R))/(max(R)-min(R));
a = [x,y];

filename = strcat('../',num2str(mu_best),num2str(G_best),'.csv');
writematrix(a, filename);
plot(x,y,'s')
%%
rmpath('src/spectral')
end