% file s_generate_compatibility.m
% brief contains script to generate golden data for IMR compatibility cases

% brief This script generates the golden data for compatibility cases
clc;
clear;
close all;

addpath('../toolchain/');
addpath('../src/forward_solver/');
addpath('../tests/');
load('file_ids.mat');

% constants / vectors
tvector     = linspace(0,12e-6,100);
threshold   = 5;
collapse    = 0;
masstrans   = 0;
vapor       = 1;
R0          = 50e-6;
Req         = R0/10;

muvec       = linspace(1e-4,1e-1,2);
Gvec        = linspace(1e2,0.25e4,2);
alphaxvec   = linspace(1e-3,1,2);
lambda1vec  = linspace(1e-7,1e-3,2);

radial_vec   = 1:4;
bubtherm_vec = 0:1;
medtherm_vec = 0:1;
stress_vec   = 0:5;

% indexing: keep only valid combinations
dims       = [numel(lambda1vec) numel(alphaxvec) numel(Gvec) numel(muvec) ...
    numel(stress_vec) numel(medtherm_vec) numel(bubtherm_vec) ...
    numel(radial_vec)];
total_full = prod(dims);

idx_all = 1:total_full;
[~,~,~,~,~,med_i,bub_i,~] = ind2sub(dims, idx_all);
valid_idx = idx_all(~(bubtherm_vec(bub_i) == 0 & medtherm_vec(med_i) == 1));
total_valid = numel(valid_idx);

filenames_fd = cell(total_valid,1);
filenames_sp = cell(total_valid,1);
for idx = 1:total_valid
    filenames_fd{idx} = sprintf('../tests/%s.mat', ids{idx});
    filenames_sp{idx} = sprintf('../tests/%s.mat', ids{idx + total_valid});
end

% parallel pool
if isempty(gcp('nocreate'))
    parpool('local',14);
end

% dispatch
futures(total_valid, 1) = parallel.FevalFuture;
for k = 1:total_valid
    idx_full = valid_idx(k);
    [lambda1_i, alphax_i, G_i, mu_i, stress_i, med_i, bub_i, rad_i] = ind2sub(dims, idx_full);
    
    P = struct( ...
        'radial',   radial_vec(rad_i), ...
        'bubtherm', bubtherm_vec(bub_i), ...
        'medtherm', medtherm_vec(med_i), ...
        'stress',   stress_vec(stress_i), ...
        'mu',       muvec(mu_i), ...
        'G',        Gvec(G_i), ...
        'alphax',   alphaxvec(alphax_i), ...
        'lambda1',  lambda1vec(lambda1_i), ...
        'tvector',  tvector, ...
        'vapor',    vapor, ...
        'collapse', collapse, ...
        'masstrans',masstrans, ...
        'Req',      Req, ...
        'R0',       R0);
    
    futures(k) = parfeval(@f_generate_goldendata_wrapper, 3, k, P);
end

% collect & save
for i = 1:total_valid
    try
        [jobID, Rf, Rs] = fetchNext(futures);
        diff = norm(abs(Rf./Rs - 1), 2);
        if norm(abs(Rf./Rs - 1), 2) < threshold
            f_savefile_fd(filenames_fd{jobID}, Rf);
            f_savefile_sp(filenames_sp{jobID}, Rs);
            fprintf('✓ Saved index %d, diff %+.5E\n', i, diff);
        else
            error('Mismatch at idx %d, diff %+.5E\n', i, diff);
        end
    catch
        fprintf('✗ Job %d failed: %s\n', i, ME.message);
        error('stopping');
    end
end

delete(gcp('nocreate'));

% helper functions
function [Rf, Rs, jobID] = f_generate_goldendata_wrapper(jobID, P)
    args = {'progdisplay',0,...
        'radial',P.radial,...
        'bubtherm',P.bubtherm, ...
        'tvector',P.tvector,...
        'vapor',P.vapor,...
        'medtherm',P.medtherm,...
        'masstrans',P.masstrans,...
        'collapse',P.collapse,...
        'lambda2',0, ...
        'Req',P.Req,...
        'R0',P.R0,...
        'mu',P.mu,...
        'G',P.G,...
        'alphax',P.alphax, ...
        'lambda1',P.lambda1,...
        'stress',P.stress};
    [~, Rf] = f_imr_fd(args{:}, 'Nt', 70, 'Mt', 70);
    [~, Rs] = f_imr_spectral(args{:}, 'Nt', 12, 'Mt', 12);
end

function f_savefile_fd(filename, data)
    Rf = data;
    save(filename, 'Rf', '-v7');
end

function f_savefile_sp(filename, data)
    Rs = data;
    save(filename, 'Rs', '-v7');
end
