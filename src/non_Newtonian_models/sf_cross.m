function f = sf_cross(a,lambda,gammadot_R)
    f = 1/(1+(lambda*gammadot_R).^a);
end