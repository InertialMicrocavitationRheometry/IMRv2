# `f_call_params.m` Documentation

## Overview
`f_call_params` sets up the simulation environment for the Inertial Microcavitation Rheometry (IMR) solver.
It loads a case file or default case, applies user overrides, performs non-dimensionalization, checks inputs, and packages all simulation parameters needed for the solver.

## Function Signature
```matlab
[eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] = f_call_params(varargin)
```

---

## Inputs

The function accepts **parameter-value pairs** through `varargin`. Key recognized parameters are:

### Case file and basic control
- `casefile`: Path to a custom case `.m` file.
- `collapse`: Enable special initial condition calculation for collapse.

### Time configuration
- `tfin`: Final simulation time.
- `tvector`: Full time vector for simulation.

### Radial model settings
- `radial`: Radial dynamics model (integer 1–4).
- `bubtherm`: Enable bubble thermal model (0 or 1).
- `medtherm`: Enable medium thermal model (0 or 1).

### Material model settings
- `stress`: Stress/constitutive model (0–7).
- `eps3`: Enable second-order strain (0 or 1).
- `vapor`: Enable vapor phase effects (0 or 1).
- `masstrans`: Enable mass transfer (0 or 1).

### Solver options
- `method`, `spectral`, `divisions`, `nv`, `nt`, `mt`, `lv`, `lt`: Solver configuration settings.

### Initial conditions
- `r0`: Initial radius.
- `u0`: Initial wall velocity.
- `req`: Equilibrium radius.
- `p0`: Initial bubble pressure.
- `stress0`: Initial stress field value.

### Output options
- `dimout`: Enable dimensionalized outputs (0 or 1).
- `progdisplay`: Display simulation progress (0 or 1).

### Acoustic/environmental properties
- `rho8`: Far-field density.
- `gam`: Grüneisen parameter.
- `nstate`: EOS exponent.
- `p8`: Far-field pressure.
- `c8`: Far-field sound speed.
- `hugoniot_s`: Hugoniot slope.

### Waveform configuration
- `wave_type`, `pa`, `omega`, `tw`, `dt`, `mn`: Driving waveform parameters.

### Stress/viscoelastic model parameters
- `mu`, `g`, `lambda1`, `lambda2`, `alphax`, `surft`: Viscosity, elasticity, and surface tension parameters.

### Non-Newtonian viscosity options
- `du`, `mu0`, `v_a`, `v_nc`, `v_lambda`, `nu_model`: Advanced viscosity settings.

### Thermal properties
- `t8`: Far-field temperature.
- `kappa`: Thermal diffusivity.
- `atg`, `btg`, `atv`, `btv`: Temperature conductivity coefficients.
- `km`: Medium thermal conductivity.
- `dm`: Medium diffusivity.

### Mass transfer parameters
- `dmass`, `lheat`, `rv`, `ra`: Mass transfer and gas constants.

### Pressure properties
- `pv`: Vapor pressure.

---

## Outputs

The function returns the following grouped outputs:

- `eqns_opts`: Flags for physical models (radial, thermal, stress, mass transfer).
- `solve_opts`: Solver settings (method, divisions, spectral options).
- `init_opts`: Initial conditions in non-dimensional form (radius, velocity, pressure).
- `init_stress`: Initial stress tensor (depends on stress model).
- `tspan_opts`: Simulation time span or vector (non-dimensional).
- `out_opts`: Output and display control options.
- `acos_opts`: Acoustic medium properties (sound speed, compressibility).
- `wave_opts`: Driving waveform parameters and optional preloaded data.
- `sigma_opts`: Viscoelastic material properties and dimensionless numbers.
- `thermal_opts`: Non-dimensional thermal transport coefficients.
- `mass_opts`: Mass transfer parameters and non-dimensional vapor concentration.

---

## Notes

- `f_call_params` must be called before running the forward solver.
- Critical inputs must be matched correctly; otherwise, the function throws errors.
- Defaults from `default_case.m` are loaded if no `casefile` is specified.
- Some helper functions (`f_pvsat.m`, `f_init_stress.m`) are required for full execution.

