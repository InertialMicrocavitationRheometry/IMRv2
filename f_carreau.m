function [mu] = f_carreau(vmaterial,gammadot)
	mu_inf = 0.00345; 
	mu_o = 0.056; 
    if strcmp('water',vmaterial)==1
        mu = 8.9*10E-4*ones(size(gammadot));
        
    elseif strcmp('mu_0', vmaterial) == 1
        mu = mu_o*isfinite(gammadot);
    elseif strcmp('mu_inf', vmaterial) == 1
        mu = mu_inf*isfinite(gammadot);
        
    elseif strcmp('blood',vmaterial)==1
        lambda = 3.313; 
        nc = 0.3568; 
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2); 
       
    elseif strcmp('lsq_blood',vmaterial)==1
        lambda = 5.607570991983291;
        nc = 0.383559338674208;
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
       
    elseif strcmp('lsq_blood_ska',vmaterial)==1
        lambda = 4.582273132970104;
        nc = 0.218632645992341;
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);

	elseif strcmp('lsq_blood_biro',vmaterial)==1
        lambda = 2.963374935764122;
        nc = 0.407152928274111;
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
        
    elseif strcmp('lsq_blood_meri',vmaterial)==1
        lambda = 9.670808761639107;
        nc = 0.205785916818776;
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);        
        
    elseif strcmp('polystyrene', vmaterial) ==1
        mu_inf = 0; 
        muo = 4*10^6;
        lambda = 46.4; 
        nc = 0.4;
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
            
    elseif strcmp('aluminum soap', vmaterial) == 1
        mu_inf = 0.01;
        muo = 89.6; 
        lambda = 1.41; 
        nc = 0.2; 
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
        
    elseif strcmp('p-oxide', vmaterial) == 1
        mu_inf = 0; 
        muo = 15.25; 
        lambda = 1.1876; 
        nc = 0.4133; 
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
        
    elseif strcmp('h-cellulose', vmaterial) == 1
        mu_inf = 0; 
        muo = 0.22;
        lambda = 0.0664; 
        nc = 0.5088; 
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
    end
end