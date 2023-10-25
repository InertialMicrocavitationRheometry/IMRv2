function sofr = f_gammadot_r(r,R,Rdot,N,M)
% shear as a function of r (radial coordinate) calculation
sofr = zeros(N,M);
    for i = 1:N
        for j = 1:M
            if r(i,j) >=  R(i,1)
                sofr(i,j) = -2*Rdot(i,1)*(R(i,1)^2)/((r(i,j))^3);
            else
                sofr(i,j) = NaN;
            end
        end
    end
end