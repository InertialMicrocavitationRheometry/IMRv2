function mu = carreau(vmaterial,gammadot)

    if strcmp('water',vmaterial)==1
        mu = 8.9*10E-4;
        
    elseif strcmp('mu_knot', vmaterial) == 1
        mu = 0.056; 
        
    elseif strcmp('mu_inf', vmaterial) == 1
        mu = 0.0345;
        
    elseif strcmp('lsq_mu_knot', vmaterial) == 1
        mu = 0.112; 
        
    elseif strcmp('lsq_mu_inf', vmaterial) == 1
        mu = 0.002679;
        
    elseif strcmp('blood',vmaterial)==1
        mu_inf = 0.0345; 
        muo = 0.056; 
        lambda = 3.313; 
        nc = 0.3568; 
       mu = mu_inf + (muo - mu_inf).*(1+(lambda).^2*(gammadot).^2).^((nc-1)./2); 
       
    elseif strcmp('lsq_blood',vmaterial)==1
        mu_inf = 0.002679;
        mu_o = 0.112;
        lambda = 40.12;
        nc = 0.5074;
        mu = mu_inf + (mu_o - mu_inf).*(1+(lambda).^2*(gammadot).^2).^((nc-1)./2);
       
    elseif strcmp('polystyrene', vmaterial) ==1
        mu_inf = 0; 
        muo = 4*10^6;
        lambda = 46.4; 
        nc = 0.4;
        
        mu = mu_inf + (muo - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
        
        
    elseif strcmp('aluminum soap', vmaterial) == 1
        mu_inf = 0.01;
        muo = 89.6; 
        lambda = 1.41; 
        nc = 0.2; 
        
        mu = mu_inf + (muo - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
        
        
    elseif strcmp('p-oxide', vmaterial) == 1
        mu_inf = 0; 
        muo = 15.25; 
        lambda = 1.1876; 
        nc = 0.4133; 
        
        mu = mu_inf + (muo - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
        
    elseif strcmp('h-cellulose', vmaterial) == 1
        mu_inf = 0; 
        muo = 0.22;
        lambda = 0.0664; 
        nc = 0.5088; 
        
        mu = mu_inf + (muo - mu_inf).*(1+(lambda).^2.*(gammadot).^2).^((nc-1)./2);
      
    end

end