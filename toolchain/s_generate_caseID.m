clc;
clear;
close all;
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789';
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 7;
%generate random string
ntests = 100000;
ids = cell(ntests,1);
for i = 1:ntests
    filename = s( ceil(rand(1,sLength)*numRands) );
    ids{i} = filename;
end
% Find the last non-empty cell index
lastNonEmptyIdx = find(~cellfun(@isempty, ids), 1, 'last');
% Truncate the array, keeping empty cells within range
ids = ids(1:lastNonEmptyIdx);
save('file_ids','ids');
