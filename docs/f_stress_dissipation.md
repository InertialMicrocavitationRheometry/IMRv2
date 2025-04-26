# `f_stress_dissipation.m` Documentation

## Overview
`f_stress_dissipation` calculates the mechanical dissipation term (`τ:∇u`) contributing to energy loss, accounting for both finite-difference and spectral solvers.

## Major Sections
1. **Radial Stretch Calculation**: Compute `Rst`, `x2`, `ix2`, `x4`.
2. **Finite-Difference Dissipation**: Branch for models 0–5.
3. **Spectral Dissipation**: Alternative calculation if `spectral == 1`.

## Inputs

| Name           | Description                            | Units |
|----------------|----------------------------------------|-------|
| `stress`       | Stress model selector (0–5)            | –     |
| `spectral`     | Solver type flag                       | –     |
| `Req`, `R`     | Equilibrium and current radius         | –     |
| `Rdot`         | Wall velocity                          | –     |
| `Ca`, `Br`     | Cauchy number and Brinkman number      | –     |
| `Re8`, `alphax`| Reynolds and quadratic exponent        | –     |
| `yT2`, `yT3`   | Collocation grid vectors               | –     |
| `iyT3`, `iyT4`, `iyT6` | Inverses of grid powers        | –     |
| `X`, `ZZT`     | State vectors and spectral matrix      | –     |
| `ivisco1`, `ivisco2` | Auxiliary stress indices         | –     |
| `fnu`, `DRe`   | Non-Newtonian friction and number      | –     |

## Outputs

| Name        | Description                           | Units |
|-------------|---------------------------------------|-------|
| `taugradu`  | Mechanical dissipation (`τ:∇u`)        | –     |

## Notes
- Spectral branch uses `ZZT*(X(ivisco1)-X(ivisco2))` for dissipation.
