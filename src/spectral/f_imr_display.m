% Command window display
function [] = f_imr_display(radial, bubtherm, masstrans, stress, spectral, eps3, Re8, De, Ca, LAM)
    if radial == 1
        eqn = 'Rayleigh Plesset equation';
    elseif radial == 2
        eqn = 'Keller-Miksis pressure';
    elseif radial == 3
        eqn = 'Keller-Miksis enthalpy';
    else
        eqn = 'Gilmore';
    end
    const = 'none';
    
    if bubtherm == 1
        if medtherm == 1, therm = 'full';
        else
            therm = 'cold-medium approximation';
        end
    else
        therm = 'polytropic approximation';
    end

    if masstrans == 1
        mass = 'mass transfer in the bubble';
    else 
        mass = 'no mass transfer in the bubble';
    end

    if stress == 1
        if Ca == Inf
            const = 'Newtonian fluid';
        else
            const = 'neo-Hookean Kelvin-Voigt';
        end
    elseif stress == 2
        if Ca ~= Inf && LAM == 0
            const = 'linear Zener';
        elseif Ca == Inf && LAM == 0
            const = 'linear Maxwell';
        elseif Ca == Inf && LAM ~= 0
            const = 'linear Jeffreys';
        else
            const = 'Kelvin-yangChurch series';
        end
    elseif stress == 3
        if Ca ~= Inf && LAM == 0
            const = 'upper-convective Zener';
        elseif Ca == Inf && LAM == 0
            const = 'upper-convective Maxwell';
        elseif Ca == Inf && LAM ~= 0
            const = 'Oldroyd-B';
        end
    elseif stress == 4 
        const = 'Phan-Thien-Tanner';
    else 
        const = ['Giesekus(' num2str(eps3) ')'];
    end

    if spectral == 1
        solut = 'spectral method';
    else
        solut = 'ODE formulation';
    end
    
    % display run settings
    disp('--- IMRV2 SETTINGS ---');
    disp(['Radial dynamics: ' eqn]);
    disp(['Medium rheology: ' const]);
    disp(['Thermal effects: ' therm]);
    disp(['Mass effects: ' mass]);
    disp(['Solution method: ' solut]);
    disp('--- Dimensionless numbers ---');
    disp(['Re8 = ' num2str(Re8,'%10.10f')]);
    disp(['De = ' num2str(De,'%10.10f')]);
    disp(['Ca = ' num2str(Ca,'%10.10f')]);
    disp(['LM = ' num2str(LAM,'%10.10f')]);
end