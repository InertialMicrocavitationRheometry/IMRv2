function f = sf_carreau(nc,lambda,gammadot_R)
    f = (1+(lambda).^2.*(gammadot_R).^2).^((nc-1)./2);
end