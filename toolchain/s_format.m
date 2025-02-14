% toolchain/auto_indent.m - Command-line MATLAB Auto-Indentation (No Editor)
folders = ["../toolchain", "../src"];
fileList = {}; % Store file paths

% Get the name of this script to avoid modifying itself
thisScriptFile = strcat('s_format', '.m'); % e.g., 'auto_indent.m'

% Collect all .m files in the specified directories
for folder = folders
    files = dir(fullfile(folder, '**', '*.m'));
    for f = files'
        filePath = fullfile(f.folder, f.name);

        % Skip this script itself
        [~, fileName, ext] = fileparts(filePath);
        if strcmp([fileName, ext], thisScriptFile)
            fprintf('Skipping self: %s\n', filePath);
            continue;
        end
        
        fileList{end+1} = filePath; % Append file paths to cell array
    end
end

% Process each file
for i = 1:length(fileList)
    filePath = fileList{i};
    
    try
        % Read file contents
        fid = fopen(filePath, 'r');
        lines = {};
        while ~feof(fid)
            line = fgetl(fid);
            lines{end+1} = strtrim(line); % Trim whitespace
        end
        fclose(fid);
        
        % Reconstruct the script with proper indentation
        indentLevel = 0;
        formattedLines = {};
        indentStep = '    '; % Define indentation as 4 spaces (adjustable)
        multiLineContinuation = false; % Track if last line ended in "..."

        for j = 1:length(lines)
            line = lines{j};

            % Check if previous line was a continuation
            if multiLineContinuation
                formattedLines{end+1} = [repmat(indentStep, 1, indentLevel + 1), line]; % Add extra indent
                multiLineContinuation = false; % Reset flag
                continue;
            end

            % Handle "end" statements first (reduce indentation before printing)
            if startsWith(line, "end")
                indentLevel = max(indentLevel - 1, 0);
            end

            % Enforce one semicolon per line
            semicolonIdx = strfind(line, ';'); % Find all semicolons
            if length(semicolonIdx) > 1
                % Keep only the first semicolon in the original line
                newLine1 = line(1:semicolonIdx(1)); 
                newLine2 = strtrim(line(semicolonIdx(1) + 1:end)); % Everything after the first semicolon
                
                % Store first part normally
                formattedLines{end+1} = [repmat(indentStep, 1, indentLevel), newLine1];
                
                % Add remaining commands as new lines with proper indentation
                if ~isempty(newLine2)
                    formattedLines{end+1} = [repmat(indentStep, 1, indentLevel), newLine2];
                end
                continue; % Skip normal line storage since we manually split it
            end

            % Check if line is "else" or "catch" and adjust indent level
            isElseOrCatch = any(startsWith(line, ["else", "catch"]));
            if isElseOrCatch
                formattedLines{end+1} = [repmat(indentStep, 1, max(indentLevel - 1, 0)), line];
            else
                formattedLines{end+1} = [repmat(indentStep, 1, indentLevel), line];
            end

            % Increase indent level **after** processing "if", "for", "while", etc.
            if any(startsWith(line, ["if", "for", "while", "switch", "function", "try"]))
                indentLevel = indentLevel + 1;
            end

            % Detect if line ends with "..." and mark continuation
            if endsWith(line, "...")
                multiLineContinuation = true;
            end
        end
        
        % Write back the formatted content
        fid = fopen(filePath, 'w');
        fprintf(fid, '%s\n', formattedLines{:});
        fclose(fid);

        fprintf('Auto-indented: %s\n', filePath);
    catch e
        fprintf('Skipping: %s (Error: %s)\n', filePath, e.message);
    end
end

fprintf('Auto-indentation complete.\n');
