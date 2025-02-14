% toolchain/lint.m - MATLAB Lint Script
files = [dir('../toolchain/**/*.m'); dir('../src/**/*.m')];

errorFound = false; % Track if any issues are detected

for i = 1:length(files)
    
    filePath = fullfile(files(i).folder, files(i).name);
    fprintf('Checking: %s\n', filePath);
    
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
    fprintf('Linting issues detected. Failing the workflow.\n');
    exit(1); % This makes GitHub Actions fail
    else
    fprintf('No linting issues found. Passing the workflow.\n');
end
