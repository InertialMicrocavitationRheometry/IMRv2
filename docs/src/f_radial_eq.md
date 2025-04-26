# `f_radial_eq.m` Documentation

## Overview
`f_radial_eq` evaluates the radial acceleration (`Rddot`) of a bubble wall using various models (Rayleigh–Plesset, Keller–Miksis, Gilmore) and equations of state.

## Major Sections
1. **Rayleigh–Plesset (Model 1)**
2. **Keller–Miksis Pressure Form (Model 2)**
3. **Keller–Miksis Enthalpy Form (Model 3)**
4. **Gilmore with Tait EoS (Model 4)**
5. **Keller–Miksis with Mie–Grüneisen EoS (Model 5)**
6. **Gilmore with Mie–Grüneisen EoS (Model 6)**
7. **EOS Helper Functions**:
   - `f_mie_gruneisen_eos_scalar`
   - `f_mie_rho_from_p_scalar`

## Inputs

| Name        | Description                                | Units |
|-------------|--------------------------------------------|-------|
| `radial`    | Model selector (1–6)                      | –     |
| `P`         | Internal bubble pressure (non-dim)        | –     |
| `Pdot`      | Time derivative of `P`                    | –     |
| `Pf8`       | External pressure (non-dim)               | –     |
| `Pf8dot`    | Time derivative of `Pf8`                  | –     |
| `iWe`       | Inverse Weber number                      | –     |
| `R`         | Bubble radius (non-dim)                   | –     |
| `Rdot`      | Bubble wall velocity (non-dim)            | –     |
| `S`         | Stress term (non-dim)                     | –     |
| `Sdot`      | Time derivative of `S`                    | –     |
| `Cstar`     | Non-dim sound speed                       | –     |
| `sam`       | Tait EOS parameter                        | –     |
| `no`        | Tait exponent                             | –     |
| `GAMa`      | Mie–Grüneisen parameter                   | –     |
| `nstate`    | EOS exponent                              | –     |
| `nog`       | (nstate–1)/2 for MG EoS                   | –     |
| `hugoniot_s`| Hugoniot slope                            | –     |
| `JdotA`     | Viscous and thermal coupling parameter    | –     |
| `ddintfnu`  | Second derivative of non-Newtonian term   | –     |
| `iDRe`      | Inverse Debye–Stokes number              | –     |

## Outputs

| Name     | Description                             | Units |
|----------|-----------------------------------------|-------|
| `Rddot`  | Bubble wall acceleration                | –     |

## Notes
- EOS helper functions perform scalar EOS and density calculations.
- Ensure physical consistency when mixing different models.
