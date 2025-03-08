% file s_lint.m
% brief contains script run a lint checker throughout the code

% brief This script runs the Matlab linter to check for inefficiencies in
% the code
clc;
clear;

repoRoot = fileparts(mfilename('fullpath')); % Get script location
% repoRoot = '';
toolchainPath = fullfile(repoRoot, '..', 'toolchain');
srcPath = fullfile(repoRoot, '..', 'src');

files = [dir(fullfile(toolchainPath, '**', '*.m'));
dir(fullfile(srcPath, '**', '*.m'))];

errorFound = false; % Track if any issues are detected

for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    if strcmp(files(i).name,'f_multirun.m')
        continue;
    else
        fprintf('Checking: %s\n', filePath);
    end
    
    issues = checkcode(filePath, '-id');
    
    if ~isempty(issues)
        fprintf('Linting issues found in: %s\n', filePath);
        for j = 1:length(issues)
            fprintf('%s (Line %d): %s\n', issues(j).id, issues(j).line, issues(j).message);
        end
        errorFound = true; % Mark that an issue was found
    end
end

% Exit MATLAB with an error code if any issues were detected
if errorFound
    fprintf('❌ Linting issues detected. Failing the workflow.\n');
    exit(1); % Fail the workflow
else
    fprintf('✅ No linting issues found. Passing the workflow.\n');
    exit(0); % Pass the workflow
end
