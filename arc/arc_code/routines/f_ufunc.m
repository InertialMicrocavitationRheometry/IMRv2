function uofr = f_ufunc(r,R,U,N,M)
uofr = zeros(N,M);
    for i = 1:N
        for j = 1:M
            if r(i,j) >  R(i,1)
                uofr(j,i) = U(i,1).*(R(i,1)./r(i,j))^.2;
            else
                uofr(j,i) = NaN;
            end
        end
    end
end