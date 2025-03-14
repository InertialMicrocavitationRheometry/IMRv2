% file s_resolution_check.m
% brief contains script to check resolution capability of IMR

% brief This script solves IMR for different spatial resolutions
clc;
clear;
close;

addpath('../../src');

% material parameter test cases
tvector = linspace(0,100E-6,500);
threshold = 5e-2;
collapse = 1;
masstrans = 1;
vapor = 1;
count = 1;
bubtherm = 1;
medtherm = 1;
R0 = 2e-04;
Req = R0/8;
mu = 10^-1;
G = 0.125*10^4;
alphax = 10^0;
lambda1 = 10^-3;
radial = 3;
stress = 3;

% define parameter vectors and limits
Mt_vec = round(logspace(1.6990,3.6021,20));

Nt = 50;
% calculate total number of combinations
total_comb = numel(Mt_vec);
reserrorvec = cell(total_comb,1);

parpool('local',8);

parfor idx = 1:total_comb
    
    % convert linear index to subindices
    [Mt_idx] = ind2sub([numel(Mt_vec)], idx);
    Mt = Mt_vec(Mt_idx);
    
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
    
    [~,Rf] = m_imr_fd(varin{:});
    % save safely
    reserrorvec{idx} = Rf;
end
save("fdresolution_Mt.mat","reserrorvec");
