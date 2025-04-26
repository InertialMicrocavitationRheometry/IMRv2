# `m_imr_fd.m` Documentation

## Overview
`m_imr_fd` is a fourth- and sixth-order accurate finite difference solver for Rayleigh–Plesset type equations, modeling bubble dynamics with optional thermal transport and viscoelastic material response.  
It uses a non-dimensionalized formulation and is capable of handling strong nonlinearities through adaptive finite difference schemes.

---

## Major Sections of the Code

### 1. Initialization
- Calls `f_call_params(varargin{:})` to set up all problem parameters.
- Unpacks solver settings, physical properties, initial conditions, and waveform settings.

### 2. Mesh and Finite Difference Matrix Setup
- Constructs spatial grids for velocity (`v`) and stress (`s`) fields.
- Builds fourth- and sixth-order finite difference derivative matrices (using custom `f_create_matrix.m`).
- Sets up boundary condition enforcement for non-trivial PDEs.

### 3. Time Integration Loop
- Advances the solution using implicit midpoint rule and Newton iteration (`ode_imr_fd.m`).
- Updates radius, velocity, and stress fields over each time step.
- Handles possible singularities or collapses with robustness checks.

### 4. Output Assembly
- Packages time histories of radius, velocity, stress fields, temperature fields, and other diagnostics.
- Returns outputs through `varargout`, depending on simulation options.

---

## Inputs

Inputs are passed via `varargin`, routed through `f_call_params`, and then unpacked locally.

| Variable          | Description |
|-------------------|-------------|
| `radial`           | Radial dynamics model (1–4) |
| `bubtherm`         | Enable bubble thermal model (0 or 1) |
| `medtherm`         | Enable medium thermal model (0 or 1) |
| `stress`           | Stress model selection (0–7) |
| `eps3`             | Include second-order strain term (0/1) |
| `masstrans`        | Enable mass transfer (0/1) |
| `method`           | Solver method (1 = FD explicit, 2 = IMR) |
| `spectral`         | Enable spectral solver (should be 0 here) |
| `divisions`, `Nv`, `Nt`, `Mt` | Spatial and angular resolution parameters |
| `Lv`, `Lt`         | Domain lengths in velocity and angle |
| `Rzero`, `Uzero`   | Initial radius and velocity (non-dimensional) |
| `Pb_star`, `Pv_star`, `P8` | Non-dimensional pressures |
| `T8`               | Far-field temperature |
| `Req`              | Equilibrium radius |
| `alphax`           | Strain-softening parameter |
| `Szero`            | Initial stress field |
| `tspan`            | Time vector (non-dimensional) |
| `tfin`             | Final simulation time |
| `dimensionalout`   | Output in physical units (0/1) |
| `progdisplay`      | Show progress during run (0/1) |
| `Cstar`, `GAMa`, `kappa`, `nstate`, `hugoniot_s` | Acoustic properties |
| `om`, `ee`, `tw`, `dt`, `mn`, `wave_type` | Acoustic waveform parameters |
| `wave_poly`, `wave_dpoly` | Waveform polynomial fits (if needed) |
| `We`, `Re8`, `DRe` | Dimensionless fluid dynamics numbers |
| `v_a`, `v_nc`, `Ca`, `LAM`, `De`, `JdotA` | Viscoelastic model parameters |
| `nu_model`, `v_lambda_star`, `zeNO`, `iDRe` | Viscosity model settings |

---

## Outputs

`m_imr_fd` uses `varargout` to dynamically return simulation results, which may include:

| Output | Description |
|--------|-------------|
| `tspan` | Time vector (non-dimensional) |
| `R`     | Bubble radius history |
| `U`     | Wall velocity history |
| `S`     | Stress field history (if active) |
| `T`     | Temperature field (if thermal transport is active) |
| `mass fields` | Vapor/mass fraction fields (if mass transfer is active) |
| `Diagnostic info` | Optional convergence/error diagnostics |

Exact structure depends on solver options selected at runtime.

---

## Notes

- **Finite difference matrices** are constructed with optional sixth-order accuracy.
- **IMR (implicit midpoint rule)** ensures stability during strong nonlinear collapse events.
- Mesh resolution must be tuned carefully for high Weissenberg numbers or sharp waveforms.
- This function expects to work with supporting scripts: `f_create_matrix.m`, `ode_imr_fd.m`, `f_call_params.m`, and waveform data files if applicable.
