clc;
clear all; close all;
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 7;
%generate random string
ntests = 1000;
for i = 1:ntests
    filename = s( ceil(rand(1,sLength)*numRands) );
    ids{i} = filename;
end
save('file_ids','ids');
