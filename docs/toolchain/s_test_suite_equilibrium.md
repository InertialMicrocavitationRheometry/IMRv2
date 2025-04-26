# `s_test_suite_equilibrium.m` Documentation

## Overview
`s_test_suite_equilibrium.m` checks that both solvers correctly maintain equilibrium conditions (no dynamics) across model variations.

## Major Sections

1. **Initialization**
   - Clears workspace.
   - Adds paths: `toolchain/` and `src/forward_solver/`.

2. **Preallocation**
   - Defines `num_tests` based on parameter loops.
   - Initializes `errors_fd`, `errors_sp`, `failed_tests`.
   - Sets `threshold = 1e-14` for high-precision check.

3. **Parameter and Initial Condition Setup**
   - Sets fixed test parameters:  
     `masstrans = 0`, `collapse = 0`,  
     `Req = 100e-6`, `R0 = Req`, `T8 = 298.15`,  
     `tfin = 1e-5`, `tvector = linspace(0, tfin, 100)`.

4. **Nested Test Loop**
   - Iterates over `radial` (1–4), `vapor` (0/1),  
     `bubtherm` (0/1), `medtherm` (0/1), `stress` (0–5).
   - Constructs `varin` parameter sets.
   - Runs `m_imr_fd` and `m_imr_spectral` with fine meshes (`Nt=150,Mt=150` FD; `Nt=12,Mt=12` spectral).
   - Computes equilibrium error: `abs(norm(R_test,2)/10 - 1)`.
   - Records failures if error exceeds `threshold`.

5. **Summary and Exit**
   - Cleans up `failed_tests`.
   - Prints pass or fail message.
   - Exits with `0` on success or `1` on failure.

## Inputs

| Name             | Description                                    |
|------------------|------------------------------------------------|
| `m_imr_fd`       | Finite difference solver function              |
| `m_imr_spectral` | Spectral solver function                       |
| Script parameters| Equilibrium test parameters (e.g., `Req`, `R0`)|

## Outputs

- **Console output**: Equilibrium error norms and summary.
- **Exit code**: `0` if all within `threshold`, `1` otherwise.

## Notes

- Equilibrium condition implies the bubble radius remains constant.
- High-precision tolerance ensures numerical consistency.
