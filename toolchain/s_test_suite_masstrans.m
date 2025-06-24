% file s_test_suite_masstrans.m
% brief contains script to test for mass transfer condition

% brief This script runs the finite difference codes to
% test that mass transfer works for the radial, thermal and stress options

clc;
clear;
close;

addpath('../toolchain/');
addpath('../src/forward_solver/');
addpath('../tests/');
load('file_ids.mat');

num_tests = 4*1*1*6;
errors_fd = zeros(num_tests,1);
failed_tests = zeros(num_tests,1);
% define threshold
threshold = 5e-5;

fprintf('Checking L2 norm errors...\n');
shift = 4*2*2*6*2*2*2*2*2 - 4*1*1*6*2*2*2*2;

% mass transfer test case
tvector = linspace(0,50E-6,100);
masstrans = 1;
vapor = 1;
collapse = 0;
R0 = 2.0e-04;
Req = 3.5e-05;
radial_vec = 1:4;
bubtherm_vec = 1;
medtherm_vec = 1;
stress_vec = 0:5;

dims = [length(stress_vec), length(medtherm_vec), ...
    length(bubtherm_vec), length(radial_vec)];
total_combinations = prod(dims);

filenames_mass = cell(total_combinations,1);
for idx = 1:total_combinations
    filenames_mass{idx} = sprintf('../tests/%s.mat', ids{idx+shift});
end

for idx_mass = 1:total_combinations
    % indices
    [stress_idx, medtherm_idx, bubtherm_idx, radial_idx] = ind2sub(dims, idx_mass);
    
    radial   = radial_vec(radial_idx);
    bubtherm = bubtherm_vec(bubtherm_idx);
    medtherm = medtherm_vec(medtherm_idx);
    stress   = stress_vec(stress_idx);
    
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
    [~,Rm_test] = f_imr_fd(varin{:},'Nt',70,'Mt',70);
    load(filenames_mass{idx_mass});
    errors_fd(idx_mass) = norm(abs(Rm_test./Rm - 1),2);
    fprintf('Test %d: L2 norm error = %.6e\n', idx_mass, errors_fd(idx_mass));
    if (errors_fd(idx_mass) > threshold)
        failed_tests(idx_mass) = idx_mass;
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
    % exit(0);
else
    % fail the workflow
    fprintf('❌ Tests FAILED at indices: %s\n', sprintf('%d ', failed_tests));
    % exit(1);
end
