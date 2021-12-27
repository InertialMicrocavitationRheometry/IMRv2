function [l_mu] = f_char_length(R,U,lambda)

    l_mu = (2.*abs(U).*(R.^2).*lambda).^(1/3);

end