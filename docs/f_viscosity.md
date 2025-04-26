# `f_viscosity.m` Documentation

## Overview
`f_viscosity` computes non-Newtonian viscosity integrals for Carreau, Cross, Powell–Eyring, and related models, providing stress and its derivatives for bubble dynamics.

## Major Sections
1. **Shear Rate Computation**: `gammadot_R`, `gammadot_num`, `dgammadot`, `ddgammadot`.
2. **Model Selection**: Branch for `nu_model` (1–7).
3. **Integral Calculations**: Perform numerical integration with `integral`.
4. **Internal Functions**: Definitions for `sf_carreau`, `sf_carreau_yasuda`, etc.

## Inputs

| Name        | Description                            | Units   |
|-------------|----------------------------------------|---------|
| `nu_model`  | Viscosity model selector (1–7)        | –       |
| `Rdot`      | Wall velocity                          | –       |
| `R`         | Bubble radius                          | –       |
| `a`, `nc`   | Model parameters (exponents)           | –       |
| `lambda`    | Time constant for viscosity model      | s       |

## Outputs

| Name     | Description                                      | Units   |
|----------|--------------------------------------------------|---------|
| `f`      | Local viscosity factor                           | –       |
| `intf`   | Stress integral for non-Newtonian term           | –       |
| `dintf`  | First derivative of `intf`                       | –       |
| `ddintf` | Second derivative of `intf`                      | –       |

## Notes
- Use high tolerance (`1e-8`) for integrals.
- Models:
  - 1: Carreau
  - 2: Carreau–Yasuda
  - 3: Powell–Eyring
  - 4: Modified Powell–Eyring
  - 5: Cross
  - 6: Simplified Cross
  - 7: Modified Cross
