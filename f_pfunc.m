function pofr = f_pfunc(r,R,U,Udot,pA,rho8,N,M)
pofr = zeros(N,M);
% equation 4.11 in Eric Johnsen thesis
    for i = 1:N
        for j = 1:M
            if r(i,j) >  R(i,1)
                pofr(j,i) = rho8*((Udot(i,1).*R(i,1)...
                    + 2*U(i,1).^2).*(R(i,1)./r(i,j))...
                    - 0.5*U(i,1).^2.*(R(i,1)./r(i,j)).^4) + pA.*101325;
            else
                pofr(j,i) = NaN;
            end
        end
    end
end