% file s_test_suite_gaussian.m
% brief contains script to test for equilibrium conditions from a Gaussian
% pressure impulse

% brief This script runs both the finite difference and spectral codes to
% test that initial equilibrium conditions are satisfied
clc;
clear;

addpath('../toolchain/');
addpath('../src/forward_solver/');

num_tests = 4*2*2*2*5;
errors_fd = zeros(num_tests,1);
errors_sp = zeros(num_tests,1);
failed_tests = zeros(size(errors_sp));
% define threshold
threshold = 1e-2;

fprintf('Checking L2 norm errors...\n');

% equilibrium test case
masstrans = 0;
collapse = 0;
Req = 10e-6;
R0 = Req;
T8 = 298.15;
tfin = 2e-3;
tvector = linspace(0,tfin,1000);
mu = 0.5;
wave_type = 1;
dt = 3e-4;
tw = 1e-4;
pA = 1e5;
% equation options
count = 1;
for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for vapor = 0:1
                for stress = 1:5
                    if bubtherm == 0 && medtherm == 1
                        count = count + 2;
                        continue;
                    end
                    varin = {'radial',radial,...
                        'bubtherm',bubtherm,...
                        'masstrans',masstrans,...
                        'tvector',tvector,...
                        'vapor',vapor,...
                        'medtherm',medtherm,...
                        'stress',stress,...
                        'collapse',collapse,...
                        'method',23,...
                        'Req',Req,...
                        'mu',mu,...
                        'R0',R0,...
                        't8',T8,...
                        'wave_type',wave_type,...
                        'dt',dt,...
                        'tw',tw,...
                        'pA',pA};
                    [~,Rf_test] = f_imr_fd(varin{:},'Nt',50,'Mt',50);
                    errors_fd(count) = abs(1-Rf_test(end));
                    fprintf('Test %d: L2 norm error = %.6e\n', count, errors_fd(count));
                    [~,Rs_test] = f_imr_spectral(varin{:},'Nt',12,'Mt',12);
                    errors_sp(count) = abs(Rs_test(1)-Rs_test(end));
                    fprintf('Test %d: L2 norm error = %.6e\n', count+1, errors_sp(count));
                    if ( errors_fd(count)  > threshold )
                        failed_tests(count) = count;
                    end
                    if ( errors_sp(count) > threshold )
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
    % Success
    exit(0);
else
    fprintf('❌ Tests FAILED at indices: %s\n', sprintf('%d ', failed_tests));
    % Fail the workflow
    exit(1);
end
