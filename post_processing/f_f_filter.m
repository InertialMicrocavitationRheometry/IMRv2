function [feps_r] = f_f_filter(f_r,N,M)
%F_MUFILTER Summary of this function goes here
%   Detailed explanation goes here
eps = 0.01;
feps_r = zeros(size(f_r));
for i = 1:N
    for j = 1:M
        if f_r(i,j) > 1-eps
            feps_r(i,j) = NaN;
        elseif f_r(i,j) < eps
            feps_r(i,j) = NaN;
        else
            feps_r(i,j) = f_r(i,j);
        end
    end
end

end