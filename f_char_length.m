function [l_mu] = f_char_length(R,Rdot,lambda)

    l_mu = (2.*abs(Rdot).*(R.^2).*lambda).^(1/3);

end
