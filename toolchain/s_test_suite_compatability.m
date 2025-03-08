% file s_test_suite_compatability.m
% brief contains script to test for compatability conditions

% brief This script runs both the finite difference and spectral codes to
% test to ensure compabitility between including thermal conditions
clc;
clear;
addpath('../src');
load('file_ids.mat');

num_tests = 4*2*2*2*6*2;
errors_fd = zeros(num_tests,1);
errors_sp = zeros(num_tests,1);
failed_tests = zeros(size(errors_sp));
count = 1;
threshold = 5e-2;

fprintf('Checking L2 norm errors...\n');

% equation options
tvector = linspace(0,15E-6,100);

collapse = 1;
masstrans = 0;
vapor = 1;
R0 = 50e-6;
Req = R0/12;

for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for stress = 0:5
                filename1 = strcat('../tests/',ids{count+0},'.mat');
                filename2 = strcat('../tests/',ids{count+1},'.mat');
                load(filename1);
                load(filename2);
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
                    'stress',stress};
                [~,Rf_test] = m_imr_fd(varin{:},'Nt',50,'Mt',50);
                [~,Rs_test] = m_imr_spectral(varin{:},'Nt',12,'Mt',12);
                errors_fd(count) = norm(Rf-Rf_test,2);
                errors_sp(count) = norm(Rs-Rs_test,2);
                fprintf('Test %d: L2 norm error = %.6e\n', count, errors_fd(count));
                fprintf('Test %d: L2 norm error = %.6e\n', count+1, errors_sp(count));
                if (errors_fd(count) > threshold)
                    failed_tests(count) = count;
                end
                if (errors_sp(count) > threshold)
                    failed_tests(count+1) = count+1;
                end
                count = count + 2;
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
    exit(0); % Success
else
    fprintf('❌ Tests FAILED at indices: %s\n', sprintf('%d ', failed_tests));
    exit(1); % Fail the workflow
end
