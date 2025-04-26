# `default_case.m` Documentation

## Overview
`default_case.m` defines the default simulation parameters for the Inertial Microcavitation Rheometry (IMR) solver.  
It sets baseline flags for physics models, numerical solver options, initial conditions, material and environmental properties.  
These defaults are used when no custom case file is specified.

## Major Sections
1. **Physics Model Flags**:  
   - `radial`, `bubtherm`, `medtherm`, `stress`, `eps3`, `vapor`, `masstrans`

2. **Solver Options**:  
   - `TFin`, `TVector`, `method`, `spectral`, `divisions`, `Nt`, `Mt`, `Nv`, `Lv`, `Lt`

3. **Initial Conditions**:  
   - `collapse`, `R0`, `U0`, `Req`

4. **Output Options**:  
   - `dimensionalout`, `progdisplay`

5. **Acoustic Properties**:  
   - `rho8`, `GAM`, `nstate`, `hugoniot_s`, `P8`, `C8`

6. **Pressure Wave Options**:  
   - `pA`, `omega`, `TW`, `DT`, `mn`, `wave_type`

7. **Stress Model Parameters**:  
   - `S`, `G`, `lambda1`, `lambda2`, `alphax`

8. **Thermal Properties**:  
   - `kappa`, `ATg`, `BTg`, `ATv`, `BTv`, `T8`, `Km`, `Cp`

9. **Mass Transfer Properties**:  
   - `D0`, `L_heat`, `Ru`, `Rv`, `Ra`

10. **Viscosity Variables**:  
   - `mu8`, `muo`, `Dmu`, `v_a`, `v_nc`, `v_lambda`, `nu_model`

11. **Derived Pressure Variables**:  
   - `Pv`, `P0`

## Variables

| Name            | Description                                                                                       | Units         | Default / Calculation        |
|-----------------|---------------------------------------------------------------------------------------------------|---------------|------------------------------|
| `radial`        | Select radial model: 1=Rayleigh–Plesset, 2=Keller–Miksis (pressure form), 3=Keller–Miksis (enthalpy), 4=Gilmore | –             | 1                            |
| `bubtherm`      | Bubble thermal model: 0=off, 1=on                                                                  | –             | 0                            |
| `medtherm`      | Medium thermal model: 0=off, 1=on                                                                  | –             | 0                            |
| `stress`        | Constitutive model: 0=none, 1=NHKV, 2=quad. KV NH, 3=Maxwell/Jeffreys/Zener, 4=PTT, 5=Giesekus      | –             | 0                            |
| `eps3`          | Second-order strain term: 0=off, (0,0.5]=on                                                        | –             | 0                            |
| `vapor`         | Vapor pressure model: 0=ignore, 1=include                                                          | –             | 0                            |
| `masstrans`     | Mass transfer: 0=off, 1=on                                                                         | –             | 0                            |
| `TFin`          | Final simulation time                                                                              | s             | 20e-6                       |
| `TVector`       | Time vector [0, TFin]                                                                              | s             | [0 TFin]                    |
| `method`        | ODE solver: 23=ode23tb, 15=ode15s, 45=ode45                                                         | –             | 23                           |
| `spectral`      | Force spectral collocation: 0=FD, 1=spectral                                                        | –             | 0                            |
| `divisions`     | Max time steps (0=automatic)                                                                       | –             | 0                            |
| `Nt`            | Number of points inside bubble (thermal PDE)                                                       | –             | 25                           |
| `Mt`            | Number of points outside bubble (thermal PDE)                                                      | –             | 25                           |
| `Nv`            | Number of points outside bubble (stress PDE)                                                       | –             | 150                          |
| `Lv`            | Grid stretching factor (velocity domain)                                                           | –             | 3                            |
| `Lt`            | Grid stretching factor (thermal domain)                                                            | –             | 2                            |
| `collapse`      | Special collapse IC: 0=off, 1=on                                                                   | –             | 0                            |
| `R0`            | Initial bubble radius                                                                              | m             | 50e-6                        |
| `U0`            | Initial wall velocity                                                                              | m/s           | 0                            |
| `Req`           | Equilibrium radius (pre-stress bubble)                                                             | m             | R0/5                         |
| `dimensionalout`| Return outputs in physical units: 0=non-dim, 1=dim                                                   | –             | 0                            |
| `progdisplay`   | Display progress: 0=off, 1=on                                                                      | –             | 0                            |
| `rho8`          | Far-field density                                                                                  | kg/m³         | 1064                         |
| `GAM`           | Mie–Grüneisen parameter (Pa)                                                                        | Pa            | 3049.13×1e5                  |
| `nstate`        | EOS exponent                                                                                       | –             | 7.15                         |
| `hugoniot_s`    | Hugoniot slope                                                                                     | –             | 10                           |
| `P8`            | Far-field pressure                                                                                 | Pa            | 101325                       |
| `C8`            | Far-field sound speed                                                                              | m/s           | 1484                         |
| `pA`            | Pressure amplitude                                                                                 | Pa            | 0                            |
| `omega`        | Driving frequency                                                                                  | rad/s         | 0                            |
| `TW`            | Gaussian pulse width                                                                               | s             | 0                            |
| `DT`            | Pulse delay                                                                                        | s             | 0                            |
| `mn`            | Waveform power exponent                                                                            | –             | 0                            |
| `wave_type`     | Waveform type: 0=impulse,1=Gaussian,2=histotripsy,3=heaviside                                       | –             | 0                            |
| `S`             | Surface tension                                                                                     | N/m           | 0.072                        |
| `G`             | Shear modulus                                                                                      | Pa            | 1e3                          |
| `lambda1`       | Relaxation time                                                                                    | s             | 1e-7                         |
| `lambda2`       | Retardation time                                                                                   | s             | 1e-8                         |
| `alphax`        | Quadratic KV exponent                                                                              | –             | 0.25                         |
| `kappa`         | Ratio of specific heats                                                                            | –             | 1.47                         |
| `ATg`           | Thermal cond. coeff. (gas, W/m·K²)                                                                 | W/m·K²        | 5.28e-5                      |
| `BTg`           | Thermal cond. coeff. (gas, W/m·K)                                                                  | W/m·K         | 1.165e-2                     |
| `ATv`           | Thermal cond. coeff. (vapor, W/m·K²)                                                               | W/m·K²        | 3.30e-5                      |
| `BTv`           | Thermal cond. coeff. (vapor, W/m·K)                                                                | W/m·K         | 1.742e-2                     |
| `T8`            | Far-field temperature                                                                              | K             | 298.15                       |
| `Km`            | Thermal conductivity (medium)                                                                     | W/m·K         | 0.55                         |
| `Cp`            | Specific heat (medium)                                                                              | J/kg·K        | 4.181e3                      |
| `D0`            | Mass diffusivity                                                                                   | m²/s          | 24.2e-6                      |
| `L_heat`        | Latent heat of vaporization                                                                        | J/kg          | 2.26476×10^6                 |
| `Ru`            | Universal gas constant                                                                             | J/mol·K       | 8.3144598                    |
| `Rv`            | Vapor gas constant                                                                                 | J/kg·K        | Ru/18.01528e-3               |
| `Ra`            | Air gas constant                                                                                   | J/kg·K        | Ru/28.966e-3                 |
| `mu8`           | Infinite shear-rate viscosity                                                                      | Pa·s          | 8.3283e-4                    |
| `muo`           | Zero-shear viscosity                                                                               | Pa·s          | 8.3283e-4                    |
| `Dmu`           | Viscosity difference                                                                               | Pa·s          | 0                            |
| `v_a`           | Viscosity power factor 1                                                                           | –             | 0                            |
| `v_nc`          | Viscosity power factor 2                                                                           | –             | 0                            |
| `v_lambda`      | Viscosity time scale                                                                               | s             | 0                            |
| `nu_model`      | Viscosity model type                                                                               | 0=Newtonian...7=modified Cross | 0             |
| `Pv`            | Vapor pressure from `f_pvsat(T8)`                                                                  | Pa            | computed                     |
| `P0`            | Initial internal bubble pressure                                                                   | Pa            | (P8 + 2*S/Req - Pv*vapor)*(Req/R0)^3 |

## Notes
- Modify any variable here to change the default simulation behavior.
- `Pv` uses an empirical fit: see `f_pvsat.m`.
- `Req` default ensures pre-stressed bubble based on `R0`.
