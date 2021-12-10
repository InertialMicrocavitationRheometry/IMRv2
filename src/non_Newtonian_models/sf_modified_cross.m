function f = sf_modified_cross(a,nc,lambda,gammadot_R)
    f = 1/((1+(lambda*gammadot_R)^nc).^a);
end