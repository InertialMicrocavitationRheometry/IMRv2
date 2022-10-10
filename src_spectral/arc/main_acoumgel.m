clear; close all; clc;
format long;

options = optimset('TolX',1E-10,'MaxIter',50);
pA = 1; 
rho8 = 1000;
p8 = 1E5;
kappa = 1.4;
c8 = 1500;
pgo0 = 2300;
%stress model
voigt = 0;
neohook = 1;
%bubble wall model
rpe = 0;
enth = 1;
gil = 0;
poly = 1;
Sval = 0.07;
muval = 0.115;
% muval = 1E-3;
Gval = 0;
alpha0 = 1; 
epsiloncrit = 0.29;
routines = '/mnt/hdd/licbrc2020/routines/';
addpath(routines);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A. non-condensible gas pressure evaluation first %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitting = 1;
fitR0 = 1;
shift = 0;
Gval = 1.13E3;
filepre = '/mnt/hdd/licbrc2020/um/timeShifted/histotripsy3_timeShifted';
run /mnt/hdd/licbrc2020/um/m_um_acowater.m
save /mnt/hdd/licbrc2020/mat/um_acogel3R0_H.mat;
shift = 0;
% fitR0 = 2;
Gval = 21.7E3;
filepre = '/mnt/hdd/licbrc2020/um/timeShifted/histotripsy10_timeShifted';
run /mnt/hdd/licbrc2020/um/m_um_acowater.m
save /mnt/hdd/licbrc2020/mat/um_acogel10R0_H.mat;
rmpath(routines);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B. Plotting temporal data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear;
format long;
set(gcf,'renderer','Painters');
load('/mnt/hdd/licbrc2020/mat/um_acogel3R0_H.mat');
fontaxis0 = 'Times';
fontaxis0 = 'LM Roman 9';
fontaxis1 = 'Times';
fontaxis1 = 'LM Roman 9';
shift = 0;
fitting = 0;
addpath(routines);
run /mnt/hdd/licbrc2020/um/m_um_acogel.m;
save /mnt/hdd/licbrc2020/mat/um_acogel3R0_H.mat;
sfdir = '/mnt/hdd/licbrc2020/sf_aco/';
% svformat = 'epsc';
svformat = 'png';
suffix = 'H';
figure(1)
ofilename = strcat(sfdir,'um_RTR0_gel3_dim',suffix);
saveas(gcf,ofilename,svformat);
figure(2)
ofilename = strcat(sfdir,'um_RTR0_gel3_ndim',suffix);
saveas(gcf,ofilename,svformat);
figure(3)
ofilename = strcat(sfdir,'um_RTR0_gel3_LE',suffix);
saveas(gcf,ofilename,svformat);
close all; clc; 
format long;
set(gcf,'renderer','Painters');
load('/mnt/hdd/licbrc2020/mat/um_acogel10R0_H.mat');
fontaxis0 = 'Times';
fontaxis0 = 'LM Roman 9';
fontaxis1 = 'Times';
fontaxis1 = 'LM Roman 9';
shift = 0;
fitting = 0;
run /mnt/hdd/licbrc2020/um/m_um_acogel.m;
save /mnt/hdd/licbrc2020/mat/um_acogel10R0_H.mat;
sfdir = '/mnt/hdd/licbrc2020/sf_aco/';
% svformat = 'epsc';
svformat = 'png';
suffix = 'H';
figure(1)
ofilename = strcat(sfdir,'um_RTR0_gel10_dim',suffix);
saveas(gcf,ofilename,svformat);
figure(2)
ofilename = strcat(sfdir,'um_RTR0_gel10_ndim',suffix);
saveas(gcf,ofilename,svformat);
figure(3)
ofilename = strcat(sfdir,'um_RTR0_gel10_LE',suffix);
saveas(gcf,ofilename,svformat);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C. Plotting histogram data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
format long;
sfdir = '/mnt/hdd/licbrc2020/sf_aco/';
% svformat = 'epsc';
svformat = 'png';
histcount = 5;
fitR0 = 1;
shift = 0;
load('/mnt/hdd/licbrc2020/mat/um_acogel3R0_H.mat');
fontaxis0 = 'LM Roman 9';
fontaxis1 = 'Times';
addpath(routines);
run m_plothisto.m;
suffix = '3H';
figure(1)
ofilename = strcat(sfdir,'histR0R0_umgel',suffix);
saveas(gcf,ofilename,svformat);
figure(2)
ofilename = strcat(sfdir,'histRmaxR0_umgel',suffix);
saveas(gcf,ofilename,svformat);
figure(3)
ofilename = strcat(sfdir,'histfitR0_umgel',suffix);
saveas(gcf,ofilename,svformat);
figure(4)
ofilename = strcat(sfdir,'histalphaR0_umgel',suffix);
saveas(gcf,ofilename,svformat);
close all;
load('/mnt/hdd/licbrc2020/mat/um_acogel10R0_H.mat');
% svformat = 'epsc';
svformat = 'png';
fontaxis0 = 'LM Roman 9';
shift = 0;
fitR0 = 1;
histcount = 5;
run m_plothisto.m;
suffix = '10H';
figure(1)
ofilename = strcat(sfdir,'histR0R0_umgel',suffix);
saveas(gcf,ofilename,svformat);
figure(2)
ofilename = strcat(sfdir,'histRmaxR0_umgel',suffix);
saveas(gcf,ofilename,svformat);
figure(3)
ofilename = strcat(sfdir,'histfitR0_umgel',suffix);
saveas(gcf,ofilename,svformat);
figure(4)
ofilename = strcat(sfdir,'histalphaR0_umgel',suffix);
saveas(gcf,ofilename,svformat);
close all;
rmpath(routines);