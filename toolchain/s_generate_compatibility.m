% file s_generate_goldendata.m
% brief contains script to generate golden data for IMR compatibility

% brief This script generates the golden data for various test suites used
% for IMR
clc;
clear;
close;

addpath('../toolchain/');
addpath('../src');
addpath('../tests/');
load('file_ids.mat');

% test case parameters
tvector = linspace(0,15E-6,100);
threshold = 5e-2;
collapse = 1;
masstrans = 0;
vapor = 1;
R0 = 50e-6;
Req = R0/12;
muvec = linspace(10^-4,10^-1,2);
Gvec = linspace(10^2,0.25*10^4,2);
alphaxvec = linspace(10^-3,10^0,2);
lambda1vec = linspace(10^-7,10^-3,2);
radial_vec = 1:4;
bubtherm_vec = 0:1;
medtherm_vec = 0:1;
stress_vec = 0:5;

dims = [numel(lambda1vec), numel(alphaxvec), numel(Gvec), numel(muvec), ...
    numel(stress_vec), numel(medtherm_vec), numel(bubtherm_vec), ...
    numel(radial_vec)];
total_comb = prod(dims);

filenames_fd = cell(total_comb,1);
filenames_sp = cell(total_comb,1);
for idx = 1:total_comb
    filenames_fd{idx} = sprintf('../tests/%s.mat', ids{idx});
    filenames_sp{idx} = sprintf('../tests/%s.mat', ids{idx + total_comb});
end

% start parallel pool
if isempty(gcp('nocreate'))
    parpool('local',16);
end

% dispatch parallel jobs
futures(total_comb) = parallel.FevalFuture;

for idx = 1:total_comb
    [lambda1idx, alphaxidx, Gidx, muidx, stress_idx, medtherm_idx, ...
        bubtherm_idx, radial_idx] = ind2sub(dims, idx);
    
    param_struct = struct( ...
        'radial', radial_vec(radial_idx), ...
        'bubtherm', bubtherm_vec(bubtherm_idx), ...
        'medtherm', medtherm_vec(medtherm_idx), ...
        'stress', stress_vec(stress_idx), ...
        'mu', muvec(muidx), ...
        'G', Gvec(Gidx), ...
        'alphax', alphaxvec(alphaxidx), ...
        'lambda1', lambda1vec(lambda1idx), ...
        'tvector', tvector, ...
        'vapor', vapor, ...
        'collapse', collapse, ...
        'masstrans', masstrans, ...
        'Req', Req, ...
        'R0', R0);
    
    futures(idx) = parfeval(@m_generate_goldendata_wrapper, 3, idx, param_struct);
end

% collect and save results
for i = 1:total_comb
    try
        [~, Rf, Rs, idx] = fetchNext(futures);
        if norm(Rs - Rf, 2) < threshold
            savefile_fd(filenames_fd{idx}, Rf);
            savefile_sp(filenames_sp{idx}, Rs);
            fprintf('✓ Saved index %d\n', idx);
        else
            error('Mismatch exceeds threshold at idx %d\n', idx);
        end
    catch ME
        fprintf('✗ Job %d failed: %s\n', i, ME.message);
    end
end

% Optionally shut down the pool
delete(gcp('nocreate'));

function [Rf, Rs, idx_out] = m_generate_goldendata_wrapper(idx, param_struct)
    % Unpack parameter struct
    radial   = param_struct.radial;
    bubtherm = param_struct.bubtherm;
    medtherm = param_struct.medtherm;
    stress   = param_struct.stress;
    mu       = param_struct.mu;
    G        = param_struct.G;
    alphax   = param_struct.alphax;
    lambda1  = param_struct.lambda1;
    tvector  = param_struct.tvector;
    vapor    = param_struct.vapor;
    collapse = param_struct.collapse;
    masstrans = param_struct.masstrans;
    Req      = param_struct.Req;
    R0       = param_struct.R0;
    
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
    
    [~, Rf] = m_imr_fd(varin{:}, 'Nt', 70, 'Mt', 70);
    [~, Rs] = m_imr_spectral(varin{:}, 'Nt', 12, 'Mt', 12);
    idx_out = idx;
end

function savefile_fd(filename,data)
    Rf = data;
    save(filename,"Rf");
end

function savefile_sp(filename,data)
    Rs = data;
    save(filename,"Rs");
end
