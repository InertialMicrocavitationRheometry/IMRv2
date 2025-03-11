% file s_resolution_check.m
% brief contains script to check resolution capability of IMR

% brief This script solves IMR for different spatial resolutions
clc;
clear;
close;

addpath('../../src');
load('resolution.mat');

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
Nt_vec = 50:25:675;

% calculate total number of combinations
total_comb = numel(Nt_vec);
reserrorvec = zeros(total_comb,1);


for idx = 1:total_comb

    % convert linear index to subindices
    [Nt_idx] = ind2sub([numel(Nt_vec)], idx);
    Nt = Nt_vec(Nt_idx);
    
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
        'Mt',Nt};

    [~,Rf] = m_imr_fd(varin{:});
    % save safely in unique filenames
    reserrorvec(idx) = norm(Rf-Rres,2);
end
save("resolution_data.mat","reserrorvec");
% figure(1)
% hold on;
% plot(Ntvec,reserrorvec,'rs');