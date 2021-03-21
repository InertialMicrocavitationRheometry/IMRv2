function S = f_visG(nc,lambda,Rdot,R)
%F_VIS Summary of this function goes here
%   Detailed explanation goes here
    S = zeros(size(R));
    for i = 1:length(R)
        S(i) = ...
            integral(@(r) carreautau(r,Rdot(i),R(i),nc,lambda),R(i),Inf,...
            'RelTol',1E-10,'AbsTol',1E-13);
    end
end

function tau = carreautau(r,Rdot,R,nc,lambda)
    tau = (1./r.^4).*(1+...
        lambda.^2.*(4.*Rdot.^2.*R.^4./r.^6)).^((nc-1)/2);    
end