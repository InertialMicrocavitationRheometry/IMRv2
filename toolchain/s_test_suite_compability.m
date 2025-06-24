% file s_test_suite_compatibility.m
% brief contains script to test for compatibility conditions

% brief This script runs both the finite difference and spectral codes to
% test to ensure compatibility between including thermal conditions
clc;
clear;
close;

addpath('../toolchain/');
addpath('../src/forward_solver/');
addpath('../tests');
load('file_ids.mat');



fprintf('Checking L2 norm errors...\n');

% equation options
tvector = linspace(0,12E-6,100);
threshold   = 1e-4;
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

errors_fd = zeros(total_valid,1);
errors_sp = zeros(total_valid,1);
failed_tests = zeros(size(errors_sp));

for k = 1:total_valid
    idx_full = valid_idx(k);
    [lambda1_i, alphax_i, G_i, mu_i, stress_i, med_i, bub_i, rad_i] = ind2sub(dims, idx_full);
    
    varin = { ...
        'radial',   radial_vec(rad_i), ...
        'bubtherm', bubtherm_vec(bub_i), ...
        'medtherm', medtherm_vec(med_i), ...
        'stress',   stress_vec(stress_i), ...
        'mu',       muvec(mu_i), ...
        'G',        Gvec(G_i), ...
        'alphax',   alphaxvec(alphax_i), ...
        'lambda1',  lambda1vec(lambda1_i), ...
        'lambda2',  0, ...
        'tvector',  tvector, ...
        'vapor',    vapor, ...
        'collapse', collapse, ...
        'masstrans',masstrans, ...
        'Req',      Req, ...
        'R0',       R0};
    
    load(filenames_fd{k});
    [~,Rf_test] = f_imr_fd(varin{:},'Nt',70,'Mt',70);
    errors_fd(k) = norm(abs(Rf./Rf_test - 1),2);
    fprintf('Finite test %d: L2 norm error = %.6e\n', k, errors_fd(k));
    if (errors_fd(k) > threshold)
        failed_tests(k) = k;
    end
    
    load(filenames_sp{k});
    [~,Rs_test] = f_imr_spectral(varin{:},'Nt',12,'Mt',12);
    errors_sp(k) = norm(abs(Rs./Rs_test - 1),2);
    fprintf('Spectral test %d: L2 norm error = %.6e\n', k, errors_sp(k));
    if (errors_sp(k) > threshold)
        failed_tests(k+1) = k+1;
    end
    
end


% find the last non-empty cell index
lastNonEmptyIdx = find(failed_tests ~= 0, 1, 'last');
% truncate the array, keeping empty cells within range
failed_tests = failed_tests(1:lastNonEmptyIdx);

% remove zeros from failed_tests
failed_tests(failed_tests == 0) = [];

if isempty(failed_tests)
    % success
    fprintf('✅ All tests PASSED.\n');
    % success
    % exit(0);
else
    % fail the workflow
    fprintf('❌ Tests FAILED at indices: %s\n', sprintf('%d ', failed_tests));
    % failed
    % exit(1);
end
