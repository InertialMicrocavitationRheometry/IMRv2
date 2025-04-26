# `m_imr_spectral.m` Documentation

## Overview
`m_imr_spectral` is a Chebyshev spectral collocation solver for Rayleigh–Plesset type equations, modeling bubble dynamics with optional thermal transport and viscoelastic material behavior.  
It uses high-accuracy spectral differentiation matrices to discretize spatial derivatives and advance the solution robustly in time.

---

## Major Sections of the Code

### 1. Initialization
- Calls `f_call_params(varargin{:})` to load all simulation parameters.
- Unpacks solver settings, material properties, and initial conditions.

### 2. Spectral Grid and Matrix Construction
- Constructs a Chebyshev collocation grid between [-1,1].
- Builds spectral differentiation matrices for first and second derivatives.
- Scales matrices to the physical domain.
- Enforces boundary conditions analytically within the spectral framework.

### 3. Time Integration Loop
- Solves the resulting system using an implicit midpoint rule and Newton iteration.
- Advances radius, wall velocity, and (optionally) stress and temperature fields at each timestep.

### 4. Output Assembly
- Packages time histories of radius, velocity, stress, temperature, and mass transfer properties.
- Returns output variables dynamically through `varargout`.

---

## Inputs

Inputs are passed through `varargin` to `f_call_params` and unpacked locally.

| Variable          | Description |
|-------------------|-------------|
| `radial`           | Radial dynamics model (1–4) |
| `bubtherm`         | Enable bubble thermal model (0 or 1) |
| `medtherm`         | Enable medium thermal model (0 or 1) |
| `stress`           | Stress model type (0–7) |
| `eps3`             | Include second-order strain term (0/1) |
| `masstrans`        | Enable mass transfer (0/1) |
| `method`           | Solver method (spectral here) |
| `spectral`         | Spectral mode flag (should be 1) |
| `divisions`, `Nv`, `Nt`, `Mt` | Spectral resolution settings |
| `Lv`, `Lt`         | Physical domain lengths |
| `Rzero`, `Uzero`   | Initial bubble radius and wall velocity (non-dimensional) |
| `Pb_star`, `Pv_star`, `P8` | Pressures |
| `T8`               | Far-field temperature |
| `Req`              | Equilibrium radius |
| `alphax`           | Strain-softening parameter |
| `Szero`            | Initial stress tensor or value |
| `tspan`            | Non-dimensional simulation times |
| `tfin`             | Final simulation time |
| `dimensionalout`   | Output in dimensional units (0/1) |
| `progdisplay`      | Display simulation progress (0/1) |
| `Cstar`, `GAMa`, `kappa`, `nstate`, `hugoniot_s` | Acoustic medium properties |
| `om`, `ee`, `tw`, `dt`, `mn`, `wave_type` | Driving waveform settings |
| `wave_poly`, `wave_dpoly` | Waveform fits for custom shapes |
| `We`, `Re8`, `DRe` | Dimensionless fluid dynamics numbers |
| `v_a`, `v_nc`, `Ca`, `LAM`, `De`, `JdotA` | Viscoelastic material settings |
| `nu_model`, `v_lambda_star`, `zeNO`, `iDRe` | Viscosity model settings |

---

## Outputs

Outputs are dynamically assigned through `varargout`, typically:

| Output | Description |
|--------|-------------|
| `tspan` | Time vector (non-dimensional) |
| `R`     | Bubble radius history |
| `U`     | Wall velocity history |
| `S`     | Stress field history (if stress model active) |
| `T`     | Temperature profile (if thermal modeling active) |
| `mass fields` | Mass diffusion or vapor concentration fields (if enabled) |
| `Diagnostics` | Optional convergence/error data |

The actual number and type of outputs depend on solver settings and model features selected at runtime.

---

## Notes

- **Spectral accuracy** requires sufficient resolution (`Nv`, `Mt`) to capture steep gradients or oscillations.
- **Boundary conditions** are enforced analytically to prevent spurious oscillations.
- **IMR (implicit midpoint rule)** provides stability even under strong collapse dynamics.
- Requires helper functions such as `f_create_spectral_matrix.m`, `ode_imr_spectral.m`, and waveform data files if applicable.
