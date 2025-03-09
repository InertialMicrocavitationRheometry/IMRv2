% file s_format.m
% brief contains script to run formatting throughout the toolchain and src
% folders

% brief This script runs formats the toolchain and src code
close; 
clear;
clc;

folders = ["../toolchain", "../src","../docs","../graphics"];
fileList = cell(100,1);
count = 1;
thisScriptFile = strcat('s_format', '.m');

for folder = folders
    files = dir(fullfile(folder, '**', '*.m'));
    for f = files'
        filePath = fullfile(f.folder, f.name);

        [~, fileName, ext] = fileparts(filePath);
        if strcmp([fileName, ext], thisScriptFile)
            fprintf('Skipping self: %s\n', filePath);
            continue;
        end
        
        fileList{count} = filePath;
        count = count + 1;
    end
end

lastNonEmptyIdx = find(~cellfun(@isempty, fileList), 1, 'last');
fileList = fileList(1:lastNonEmptyIdx);

for i = 1:length(fileList)
    filePath = fileList{i};
    
    try
        fid = fopen(filePath, 'r');
        lines = cell(5000,1);
        count = 1;
        while ~feof(fid)
            line = fgetl(fid);
            lines{count} = strtrim(line);
            count = count + 1;
        end
        fclose(fid);

        lastNonEmptyIdx = find(~cellfun(@isempty, lines), 1, 'last');
        lines = lines(1:lastNonEmptyIdx);

        indentLevel = 0;
        formattedLines = cell(5000,1);
        count = 1;
        indentStep = '    ';
        multiLineContinuation = false;

        for j = 1:length(lines)
            line = lines{j};

            if startsWith(line, "end")
                indentLevel = max(indentLevel - 1, 0);
            end

            semicolonIdx = strfind(line, ';');
            if length(semicolonIdx) > 1
                newLine1 = line(1:semicolonIdx(1)); 
                newLine2 = strtrim(line(semicolonIdx(1) + 1:end));
                formattedLines{count} = [repmat(indentStep, 1, indentLevel), newLine1];
                count = count + 1;
                if ~isempty(newLine2)
                    formattedLines{count} = [repmat(indentStep, 1, indentLevel), newLine2];
                    count = count + 1;
                end
                continue;
            end

            if multiLineContinuation
                formattedLines{count} = [repmat(indentStep, 1, indentLevel + 1), line];
            else
                if any(startsWith(line, ["else", "catch"]))
                    formattedLines{count} = [repmat(indentStep, 1, max(indentLevel - 1, 0)), line];
                else
                    formattedLines{count} = [repmat(indentStep, 1, indentLevel), line];
                end
            end
            count = count + 1;

            multiLineContinuation = endsWith(line, "...");

            if ~multiLineContinuation && any(startsWith(line, ["if", "parfor", "for", "while", "switch", "function", "try"]))
                indentLevel = indentLevel + 1;
            end
        end

        lastNonEmptyIdx = find(~cellfun(@isempty, formattedLines), 1, 'last');
        formattedLines = formattedLines(1:lastNonEmptyIdx);

        fid = fopen(filePath, 'w');
        fprintf(fid, '%s\n', formattedLines{:});
        fclose(fid);

        fprintf('Auto-indented: %s\n', filePath);
    catch e
        fprintf('Skipping: %s (Error: %s)\n', filePath, e.message);
    end
end

fprintf('Auto-indentation complete.\n');