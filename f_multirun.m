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
%neoHook = 1; linear elasticity = 0
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
%% This is for 10% PA/0.06% BIS experimental data
addpath("../3dasm_data")
load("exp_data/PA_10%_0.06%_completed.mat");
%extract quantities of interest only for the first file
a = []; b=[];
for kk = 1:37
    file = expts(kk);
    t = file.t_norm;
    R = file.R_norm;
    %x = 2*(t./(t(end)-t(1)))- (t(end)+t(1))/(t(end)-t(1));
    %y = 2*(R./(max(R)-min(R)))-(max(R)+min(R))/(max(R)-min(R));
    %a = [a;x;y];
    tidx = sum(t>0)-1;
    Nidx = length(t);
    range = Nidx-tidx:((Nidx-tidx)+181);
    t = t(range);
    R = R(range);
    if (t(end) > 4.9 && t(end) < 5.1)
        a = [a;t;R];
    end
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

b = b';
%t0_all = b(:,1:2:end-1);
%t0_avg = mean(t0_all,2);
R0_all = b(:,2:2:end);
R0_avg = mean(R0_all,2);

%filename = strcat('file_,num2str(kk),'.csv');
%writematrix(c,"Experimental_Data_10PA_.06BIS.csv")
writematrix(c,"normalized_unscaled_Experimental_Data_10PA_.06BIS.csv")
R_std = std(R_all,1,2);

figure(1)
hold on
patch([t_avg; flip(t_avg)], [R_avg + R_std; flipud(R_avg-R_std)],'c','EdgeColor','b')
plot(t_avg,R_avg,'k')
hold off
%saveas(gcf,'exp_data_std_cloud.png')
saveas(gcf,'normized_unscaled_exp_data_std_cloud.png')
%% This is for one simulated data with the best fit G and mu
G_best = 1.39e+04;
mu_best = 0.1070;
%in f_call_params, changed TFin to 2E-4 in accordance to observations from
%experimental end times 

[t,R,~,~]=f_imrv2('linelas',0,'neoHook',1,'mu',mu_best,'g',G_best,'R0',R0_avg);
t = t./t(end);
R = R./R(1);
%x = 2*(t./(t(end)-t(1)))- (t(end)+t(1))/(t(end)-t(1));
%y = 2*(R./(max(R)-min(R)))-(max(R)+min(R))/(max(R)-min(R));
%a = [x,y];
a = [t,R];

filename = strcat('../../Simulation_Data_',num2str(mu_best),'_',num2str(G_best),'.csv');
writematrix(a, filename);
figure(2)
plot(x,y,'s')
%%
rmpath('../3dasm_data/')
rmpath('src/spectral')
end