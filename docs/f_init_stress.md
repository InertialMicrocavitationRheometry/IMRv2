# `f_init_stress.m` Documentation

## Overview
`f_init_stress` determines the initial stress by numerically solving bubble collapse dynamics with a shooting method, ensuring the stress at maximum radius matches the selected viscoelastic model.

## Major Sections
1. **Shooting Setup**: Configure `fminsearch` options (`sopts`).
2. **Initial Guess via Shooting**: Use `f_maxwell_growth_iter` to find wall velocity that maximizes stress integral.
3. **Time Integration**: Run `ode23tb` on `f_RP_bub_collapse_comp` over nondimensional time span.
4. **Stress Extraction**: Identify maximum stress from solution matrix.

### Internal Functions
- `f_maxwell_growth_iter`: Computes max radius error for shooting.
- `f_RP_bub_collapse_comp`: Nondimensional Rayleigh–Plesset ODE for collapse with stress evolution.

## Inputs

| Name      | Description                            | Units   |
|-----------|----------------------------------------|---------|
| `Req`     | Equilibrium radius                     | –       |
| `Re`      | Reynolds number                        | –       |
| `Ca`      | Cauchy number                          | –       |
| `De`      | Deborah number                         | –       |
| `We`      | Weber number                           | –       |
| `CL`      | Speed of sound ratio (`C/C_ref`)       | –       |
| `Pv_star` | Non-dimensional vapor pressure         | –       |

## Outputs

| Name | Description                      | Units   |
|------|----------------------------------|---------|
| `S0` | Initial stress at bubble maximum | –       |

## Notes
- Uses `ode23tb` for stiff ODE integration.
- Requires consistency between stress model and dimensional parameters.
