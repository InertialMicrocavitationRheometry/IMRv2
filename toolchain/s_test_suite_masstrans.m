% file s_test_suite_masstrans.m
% brief contains script to test for mass transfer condition

% brief This script runs the finite difference codes to
% test that mass transfer works for the radial, thermal and stress options

clc;
clear;

addpath('../src/forward_solver/');
addpath('../tests/');
load('file_ids.mat');

num_tests = 4*2*2*6;
errors_fd = zeros(num_tests,1);
failed_tests = zeros(num_tests,1);
% define threshold
threshold = 1e-5;

fprintf('Checking L2 norm errors...\n');
count = 1;
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

total_combinations = 4*1*1*6;

for idx_mass = 1:total_combinations
    % indices
    [stress_idx, medtherm_idx, bubtherm_idx, radial_idx] = ...
        ind2sub([length(stress_vec), length(medtherm_vec), ...
        length(bubtherm_vec), length(radial_vec)], idx_mass);
    
    radial = radial_vec(radial_idx);
    bubtherm = bubtherm_vec(bubtherm_idx);
    medtherm = medtherm_vec(medtherm_idx);
    stress = stress_vec(stress_idx);
    
    filename1 = strcat('../tests/',ids{count+shift},'.mat');
    load(filename1);
    varin = {'radial',radial,...
        'bubtherm',bubtherm,...
        'tvector',tvector,...
        'vapor',vapor,...
        'medtherm',medtherm,...
        'stress',stress,...
        'collapse',collapse,...
        'r0',R0,...
        'req',Req,...
        'masstrans',masstrans};
    [~,Rm_test] = m_imr_fd(varin{:},'Nt',30,'Mt',70);
    errors_fd(count) = norm(abs(Rm./Rm_test - 1),2);
    fprintf('Test %d: L2 norm error = %.6e\n', count, errors_fd(count));
    if (errors_fd(count) > threshold)
        failed_tests(count) = count;
    end
    count = count + 1;
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
    % error('Tests failed');
end
