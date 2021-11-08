function [munon,mugrid,sr] = f_mufilter(munon,mugrid,sr,N,M,mu_inf)
%F_MUFILTER Summary of this function goes here
%   Detailed explanation goes here
% eps = 0.00001;
eps = 0.01;
val = mu_inf;
for i = 1:N
    for j = 1:M
        if munon(i,j) > 1-eps
            munon(i,j) = NaN;
            mugrid(i,j) = NaN;
            sr(i,j) = NaN;
        elseif munon(i,j) < eps
            munon(i,j) = NaN;
            mugrid(i,j) = NaN;
            sr(i,j) = NaN;
        else
            munon(i,j) = 0;
            mugrid(i,j) = val;
        end
    end
end

end