function sofr = sfunc(r,R,Rdot,N,M)
sofr = zeros(N,M);
    for j = 1:M
        for i = 1:N
            if r(i,j) >=  R(i,1)
                sofr(i,j) = -2*Rdot(i,1)*(R(i,1)^2)/((r(i,j))^3);
            else
                sofr(i,j) = NaN;
            end
        end
    end
end