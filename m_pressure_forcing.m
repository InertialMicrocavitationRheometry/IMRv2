function [p,dp] = m_pressure_forcing(t_patt, p_patt, t_range)

    p = makima(t_patt,p_patt,t_range);
    poly_prime = fnder(poly_ML_makima,1);
    dp = ppval(poly_prime, t_range);


end