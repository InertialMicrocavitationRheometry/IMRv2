% file s_test_suite_masstrans.m
% brief contains script to test for mass transfer condition

% brief This script runs the finite difference codes to
% test that mass transfer works for the radial, thermal and stress options

clc;
clear;

addpath('../src');
load('file_ids.mat');

num_tests = 4*2*2*6*2;
errors_fd = zeros(num_tests,1);
errors_sp = zeros(num_tests,1);
failed_tests = zeros(size(errors_sp));
% define threshold
threshold = 1e-6;

fprintf('Checking L2 norm errors...\n');
count = 1;
shift = 4*2*2*6*2*2*2*2*2;

% mass transfer test case
tvector = linspace(0,50E-6,100);
masstrans = 1;
vapor = 1;
collapse = 1;
R0 = 2e-04;
Req = 3.5e-05;
for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for stress = 0:5
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
                [~,Rf_test] = m_imr_fd(varin{:},'Nt',50,'Mt',50);
                errors_fd(count) = norm(Rf-Rf_test,2);
                fprintf('Test %d: L2 norm error = %.6e\n', count, errors_fd(count));
                if (errors_fd(count) > threshold)
                    failed_tests(count) = count;
                end
                count = count + 1;
            end
        end
    end
end

% Find the last non-empty cell index
lastNonEmptyIdx = find(failed_tests ~= 0, 1, 'last');
% Truncate the array, keeping empty cells within range
failed_tests = failed_tests(1:lastNonEmptyIdx);

% Remove zeros from failed_tests
failed_tests(failed_tests == 0) = [];

if isempty(failed_tests)
    fprintf('✅ All tests PASSED.\n');
    % exit(0); % Success
else
    fprintf('❌ Tests FAILED at indices: %s\n', sprintf('%d ', failed_tests));
    % exit(1); % Fail the workflow
end
