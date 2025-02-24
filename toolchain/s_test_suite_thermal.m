% MATLAB script to check L2 norm errors after unit tests
clc;
clear;
addpath('../src');
load('file_ids.mat');

num_tests = 4*2*2*2*6*2;
errors_fd = zeros(num_tests,1);
errors_sp = zeros(num_tests,1);
failed_tests = zeros(size(errors_sp));
threshold = 1E-2; % Define threshold

fprintf('Checking L2 norm errors...\n');

% equation options
tvector = linspace(0,15E-6,200);
count = 1;
for radial = 1:4
    for vapor = 0:1
        for bubtherm = 0:1
            for medtherm = 0:1
                for stress = 0:5
                    filename1 = strcat('../tests/',ids{count+0},'.mat');
                    filename2 = strcat('../tests/',ids{count+1},'.mat');
                    load(filename1);
                    load(filename2);
                    varin = {'radial',radial,...
                        'bubtherm',bubtherm,...
                    'tvector',tvector,...
                        'vapor',vapor,...
                    'medtherm',medtherm,...
                        'stress',stress};
                    [~,Rf_test] = m_imrv2_finitediff(varin{:},'Nt',100,'Mt',100);
                    [~,Rs_test] = m_imrv2_spectral(varin{:},'Nt',12,'Mt',12);
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
