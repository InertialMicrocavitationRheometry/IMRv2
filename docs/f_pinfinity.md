# `f_pinfinity.m` Documentation

## Overview
`f_pinfinity` calculates the time-dependent external pressure (`p8`) and its derivative (`p8dot`) that drive bubble dynamics, supporting multiple waveform types.

## Major Sections
1. **Input Unpacking**: Extracts waveform parameters from `vararg`.
2. **Custom Waveform Handling**: If `wave_type < 0`, uses polynomial fits (`wave_poly`, `wave_dpoly`).
3. **Predefined Waveforms**:
   - **Impulse** (`wave_type = 0`)
   - **Gaussian** (`wave_type = 1`)
   - **Histotripsy** (`wave_type = 2`)
   - **Heaviside Impulse** (`wave_type = 3`)
4. **Inner Functions**: Define `impulse`, `gaussian`, `histo`, and `heaviside_impulse`.

## Inputs

| Name        | Description                                     | Units   |
|-------------|-------------------------------------------------|---------|
| `t`         | Current time                                   | s       |
| `vararg`    | Vector `[om, ee, tw, dt, mn, wave_type]`       | –       |

## Outputs

| Name     | Description                                 | Units   |
|----------|---------------------------------------------|---------|
| `p8`     | External pressure at time `t`               | Pa      |
| `p8dot`  | Time derivative of external pressure        | Pa/s    |

## Notes
- `wave_poly` and `wave_dpoly` must be defined in the caller’s workspace for custom waveforms.
- Gaussian and histotripsy waveforms depend on `ee` (amplitude), `om` (frequency), `tw` (width), `dt` (delay), and `mn` (exponent).
