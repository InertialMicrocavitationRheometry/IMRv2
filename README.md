## Welcome to the Inertial Microcavitation Rheometry Version 2

### Developer(s): 
* Mauro Rodriguez, mauro_rodriguez@brown.edu
* Mirelys Carcana-Barbosa, mirelys_carcana-barbosa@brown.edu
* Original Authors: Matthew Warnez & Carlos Barajas


## Description: 
This code is a reduced version of the IMR code taken from Estrada et al. (2018) JMPS. Additional physics have been added including the Keller-Miksis with enthalpy and non-Newtonian viscosity.

## Inputs:
* tspan - time to run simulation
* R0 - Initial Radii
* NT - number of nodes for temperature and concentration fields
* NTM - number of nodes for temperature in the medium
* Pext_type - type of external pressure ('sn' = sine, 'RC' = Rayleigh collapse, 'RG' = Rayleigh growth, impulse 'ip', non-equlibrium initial conditions (laser caviation and Flynn(1975) ) 'IC'
* Pext_Amp_Freq - amplitude and frequency of external pressure [amp w]

Note: For this code the out-of-equilibrium Rayleigh Collapse the intial mass in the bubble and radii are specified

## FOR THE FOLLOWING INPUTS 0 = FALSE AND 1 = TRUE
* disptime - Displays elapsed time on the command window
* Tgrad - Models temperature gradients in the buuble
* Tmgrad- Models temperature gradients outside the buuble
* Cgrad - Models concentration gradients in the buuble

## Outputs:
* t - time vector
* T_Bubble - Temperature inside the bubble
* T_Medium - Temperature outside the bubble
* R - Bubble Radius
* U - Bubble velocity
* P - Internal bubble pressure
* C - Vapor Concentration in the bubble
* Tm - Temperature in the medium
* Dim - outputs variables in dimensional form
* Comp - 0 (ignores compressibility effects) or 1 (uses Keller- Miksis)
* Reduced - 0 utilizes full model or 1 uses Preston's reduced order model
