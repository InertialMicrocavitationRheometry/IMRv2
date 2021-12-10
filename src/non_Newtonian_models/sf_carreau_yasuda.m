function f = sf_carreau_yasuda(a,nc,lambda,gammadot_R)
    f = (1+(lambda).^a.*(gammadot_R).^a).^((nc-1)./a);
end