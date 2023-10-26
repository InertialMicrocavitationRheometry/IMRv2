function f = sf_powell_eyring(nc,lambda,gammadot_R)
    f = sinh(lambda*gammadot_R)/(lambda*gammadot_R)^nc;
end