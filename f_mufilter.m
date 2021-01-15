function [munon,mugrid,sr] = f_mufilter(munon,mugrid,sr,N,M)
%F_MUFILTER Summary of this function goes here
%   Detailed explanation goes here

for i = 1:N
    for j = 1:M
        if munon(i,j) > 0.9
            munon(i,j) = NaN;
            mugrid(i,j) = NaN;
            sr(i,j) = NaN;
        elseif munon(i,j) < 0.1
            munon(i,j) = NaN;
            mugrid(i,j) = NaN;
            sr(i,j) = NaN;
        end
    end
end

end