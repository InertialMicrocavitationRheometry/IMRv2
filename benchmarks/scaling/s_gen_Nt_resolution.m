% file s_resolution_check.m
% brief contains script to check resolution capability of IMR

% brief This script solves IMR for different spatial resolutions
clc;
clear;
close;
format long;

addpath('../../src/forward_solver/');

% material parameter test cases
tvector = linspace(0,350E-6,500);
collapse = 1;
vapor = 1;
bubtherm = 1;
medtherm = 1;
masstrans = 1;
R0 = 200e-6;
Req = R0/10;
mu = 1E-3;
G = 1E3;
alphax = 0.5;
lambda1 = 1E-3;
lambda2 = 0;
radial = 3;
stress = 2;
maxNumCompThreads(16);
% define parameter vectors and limits
Mt = 1024;
Nt_vec = 2.^(2:9);
total_comb = numel(Nt_vec);

% setup output containers
tvec = cell(total_comb,1);
Rvec = cell(total_comb,1);

% launch asynchronous jobs
for idx = 1:total_comb
    Nt = Nt_vec(idx);
    varin = {'progdisplay',0,...
        'radial',radial,...
        'bubtherm',bubtherm,...
        'tvector',tvector,...
        'vapor',vapor,...
        'medtherm',medtherm,...
        'masstrans',masstrans,...
        'collapse',collapse,...
        'lambda2',0,...
        'Req',Req,...
        'R0',R0,...
        'mu',mu,...
        'G',G,...
        'alphax',alphax,...
        'lambda1',lambda1,...
        'stress',stress,...
        'Nt',Nt,...
        'Mt',Mt};
    [tvec, Rvec] = f_imr_fd(varin{:});
    filename = strcat('fd_Nt_',num2str(Nt),'_Mt_',num2str(Mt));
    save(filename,"tvec","Rvec");
end
