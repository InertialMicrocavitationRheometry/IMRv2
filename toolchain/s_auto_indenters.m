% toolchain/auto_indent.m - Command-line MATLAB Auto-Indentation
folders = ["../toolchain", "../src"];
fileList = {}; % Store file paths

% Collect all .m files in the specified directories
for folder = folders
    files = dir(fullfile(folder, '**', '*.m'));
    for f = files'
        fileList{end+1} = fullfile(f.folder, f.name); % Append file paths to cell array
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
        
        % Reconstruct the script with proper indentation (Basic approach)
        indentLevel = 0;
        formattedLines = {};
        indentStep = '    '; % Define indentation as 4 spaces (adjustable)
        
        for j = 1:length(lines)
            line = lines{j};
            
            % Reduce indent level if line is an "end" statement
            if startsWith(line, 'end')
                indentLevel = max(indentLevel - 1, 0);
            end
            
            % Apply indentation
            formattedLines{end+1} = [repmat(indentStep, 1, indentLevel), line];
            
            % Increase indent level if the line starts a block
            if any(startsWith(line, ["if", "for", "while", "switch", "function", "try"]))
                indentLevel = indentLevel + 1;
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
