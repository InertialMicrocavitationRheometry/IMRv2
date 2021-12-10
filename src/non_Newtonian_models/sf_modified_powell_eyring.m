function f = sf_modified_powell_eyring(nc,lambda,gammadot_R)
    f = log(lambda*gammadot_R+1)/(lambda*gammadot_R)^nc;
end