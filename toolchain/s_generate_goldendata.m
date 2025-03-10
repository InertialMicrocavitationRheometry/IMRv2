% file s_generate_goldendata.m
% brief contains script to generate golden data for IMR

% brief This script generates the golden data for various test suites used
% for IMR
clc;
clear;
close;

addpath('../toolchain/');
addpath('../src');
addpath('../tests/');
load('file_ids.mat');

% material parameter test cases
tvector = linspace(0,15E-6,100);
threshold = 5e-2;
collapse = 1;
masstrans = 0;
vapor = 1;
count = 1;
R0 = 50e-6;
Req = R0/12;
muvec = linspace(10^-4,10^-1,2);
Gvec = linspace(10^2,0.25*10^4,2);
alphaxvec = linspace(10^-3,10^0,2);
lambda1vec = linspace(10^-7,10^-3,2);

% define parameter vectors and limits
radial_vec = 1:4;
bubtherm_vec = 0:1;
medtherm_vec = 0:1;
stress_vec = 0:5;
mu_len = length(muvec);
G_len = length(Gvec);
alphax_len = length(alphaxvec);
lambda1_len = length(lambda1vec);

% calculate total number of combinations
total_comb = numel(radial_vec)*numel(bubtherm_vec)*numel(medtherm_vec)*...
    numel(stress_vec)*mu_len*G_len*alphax_len*lambda1_len;

filenames_fd = cell(total_comb,1);
filenames_sp = cell(total_comb,1);
for idx = 1:total_comb
    filenames_fd{idx} = sprintf('../tests/%s.mat', ids{idx});
    filenames_sp{idx} = sprintf('../tests/%s.mat', ids{idx+total_comb});
end

parpool('local',4);

parfor idx = 1:total_comb
    % convert linear index to subindices
    [lambda1idx, alphaxidx, Gidx, muidx, stress_idx, medtherm_idx, ...
        bubtherm_idx, radial_idx] = ...
        ind2sub([lambda1_len, alphax_len, G_len, mu_len, numel(stress_vec),...
        numel(medtherm_vec), numel(bubtherm_vec), numel(radial_vec)], idx);
    
    radial = radial_vec(radial_idx);
    bubtherm = bubtherm_vec(bubtherm_idx);
    medtherm = medtherm_vec(medtherm_idx);
    stress = stress_vec(stress_idx);
    mu = muvec(muidx);
    G = Gvec(Gidx);
    alphax = alphaxvec(alphaxidx);
    lambda1 = lambda1vec(lambda1idx);
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
        'stress',stress};
    [tf,Rf] = m_imr_fd(varin{:},'Nt',70,'Mt',70);
    [ts,Rs] = m_imr_spectral(varin{:},'Nt',12,'Mt',12);
    % save safely in unique filenames
    if (norm(Rs-Rf,2) < threshold)
        disp('----> SUCCESS! <------');
        savefile_fd(filenames_fd{idx},Rf);
        savefile_sp(filenames_sp{idx},Rs);
    else
        % figure(count)
        % hold on;
        % plot(ts,Rs,'r--s');
        % plot(tf,Rf,'k-.^');
        error('error radial not working');
    end
end

% mass transfer test case
tvector = linspace(0,50E-6,100);
masstrans = 1;
vapor = 1;
collapse = 1;
R0 = 2e-04;
Req = 3.5e-05;
radial_vec = 1:4;
bubtherm_vec = 0:1;
medtherm_vec = 0:1;
stress_vec = 0:5;

% precompute combinations
total_combinations = length(radial_vec)*length(bubtherm_vec)*...
    length(medtherm_vec)*length(stress_vec);

filenames = cell(total_combinations,1);
for idx = 1:total_combinations
    filenames{idx} = sprintf('../tests/%s.mat', ids{idx+2*total_comb});
end

parfor idx = 1:total_combinations
    % indices
    [stress_idx, medtherm_idx, bubtherm_idx, radial_idx] = ...
        ind2sub([length(stress_vec), length(medtherm_vec), ...
        length(bubtherm_vec), length(radial_vec)], idx);
    
    radial = radial_vec(radial_idx);
    bubtherm = bubtherm_vec(bubtherm_idx);
    medtherm = medtherm_vec(medtherm_idx);
    stress = stress_vec(stress_idx);
    
    varin = {'radial', radial, ...
        'bubtherm', bubtherm, ...
        'tvector', tvector, ...
        'vapor', vapor, ...
        'medtherm', medtherm, ...
        'stress', stress, ...
        'collapse', collapse, ...
        'r0', R0, ...
        'req', Req, ...
        'masstrans', masstrans};
    [~, Rf] = m_imr_fd(varin{:}, 'Nt', 50, 'Mt', 50);
    % save safely in unique filenames
    savefile_fd(filenames{idx},Rf);
end

function savefile_fd(filename,data)
    Rf = data;
    save(filename,"Rf");
end

function savefile_sp(filename,data)
    Rs = data;
    save(filename,"Rs");
end
