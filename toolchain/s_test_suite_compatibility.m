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

num_tests = 4*2*2*6*2*2*2*2 - 4*1*1*6*2*2*2*2;
errors_fd = zeros(num_tests,1);
errors_sp = zeros(num_tests,1);
failed_tests = zeros(size(errors_sp));

fprintf('Checking L2 norm errors...\n');

% equation options
tvector = linspace(0,12E-6,100);
threshold = 1e-4;
collapse = 0;
masstrans = 0;
vapor = 1;
count = 1;
R0 = 50e-6;
Req = R0/10;
muvec = linspace(10^-4,10^-1,2);
Gvec = linspace(10^2,0.25*10^4,2);
alphaxvec = linspace(10^-3,10^0,2);
lambda1vec = linspace(10^-7,10^-3,2);

for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for stress = 0:5
                for muidx = 1:2
                    for Gidx = 1:2
                        for alphaxidx = 1:2
                            for lambda1idx = 1:2
                                if bubtherm == 0 && medtherm == 1
                                    continue;
                                end
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
                                
                                filename1 = strcat('../tests/',ids{count+0},'.mat');
                                load(filename1);
                                [~,Rf_test] = f_imr_fd(varin{:},'Nt',70,'Mt',70);
                                errors_fd(count) = norm(abs(Rf./Rf_test - 1),2);
                                fprintf('Finite test %d: L2 norm error = %.6e\n', count, errors_fd(count));
                                if (errors_fd(count) > threshold)
                                    failed_tests(count) = count;
                                end
                                
                                filename2 = strcat('../tests/',ids{count+num_tests},'.mat');
                                load(filename2);
                                [~,Rs_test] = f_imr_spectral(varin{:},'Nt',12,'Mt',12);
                                errors_sp(count) = norm(abs(Rs./Rs_test - 1),2);
                                fprintf('Spectral test %d: L2 norm error = %.6e\n', count, errors_sp(count));
                                if (errors_sp(count) > threshold)
                                    failed_tests(count+1) = count+1;
                                end
                                
                                count = count + 1;
                            end
                        end
                    end
                end
            end
        end
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
    exit(0);
else
    % fail the workflow
    fprintf('❌ Tests FAILED at indices: %s\n', sprintf('%d ', failed_tests));
    error('Test failed');
    % failed
    exit(1);
end
