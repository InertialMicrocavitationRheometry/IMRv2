function S = f_vis(mu_inf,mu_o,nc,lambda,Rdot,R)
%F_VIS Summary of this function goes here
%   Detailed explanation goes here
    S = zeros(size(R));
    for i = 1:length(R)
        S(i) = 4*mu_inf*Rdot(i)./R(i)+...
            integral(@(r) carreautau(r,Rdot(i),R(i),mu_inf,mu_o,nc,lambda),...
            R(i),Inf,'RelTol',1E-8,'AbsTol',1E-12);
    end
end

function tau = carreautau(r,Rdot,R,mu_inf,mu_o,nc,lambda)
    tau = (12*(Rdot.*R.^2)./r.^4).*(mu_o-mu_inf).*(1+...
        lambda.^2.*(4.*Rdot.^2.*R.^4./r.^6)).^((nc-1)/2);    
end