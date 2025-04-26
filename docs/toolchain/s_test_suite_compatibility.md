# `s_test_suite_compatibility.m` Documentation

## Overview
`s_test_suite_compatibility.m` verifies that the finite difference (`m_imr_fd`) and spectral (`m_imr_spectral`) solvers produce compatible results under various model configurations, including thermal effects.

## Major Sections

1. **Initialization**
   - Clears workspace (`clc`, `clear`, `close`).
   - Adds necessary paths: `toolchain/`, `src/forward_solver/`, and `tests/`.
   - Loads test identifiers from `file_ids.mat`.

2. **Preallocation**
   - Computes `num_tests` as the product of parameter combinations.
   - Initializes `errors_fd`, `errors_sp`, and `failed_tests` arrays.

3. **Parameter Setup**
   - Defines parameter vectors:  
     - `muvec`: viscosity values  
     - `Gvec`: shear modulus values  
     - `alphaxvec`: strain-softening exponents  
     - `lambda1vec`: relaxation times  
   - Sets error threshold (`threshold = 1e-4`) and counter (`count = 1`).

4. **Nested Test Loop**
   - Loops over combinations of:  
     `radial` (1–4), `bubtherm` (0/1), `medtherm` (0/1), `stress` (0–5),  
     `mu`, `G`, `alphax`, `lambda1` indices.
   - For each combination:
     1. Constructs `varin` parameter-value pairs.
     2. Loads reference results: `Rf` and `Rs` from test MAT files.
     3. Runs `m_imr_fd` and compares output `Rf_test`; computes L2 error.
     4. Runs `m_imr_spectral` and compares output `Rs_test`; computes L2 error.
     5. Records failures where error exceeds `threshold`.

5. **Summary and Exit**
   - Truncates and cleans `failed_tests`.
   - If all tests passed, prints success and calls `exit(0)`.
   - Otherwise, prints failed test indices and raises an error.

## Inputs

| Name             | Description                                    |
|------------------|------------------------------------------------|
| `file_ids.mat`   | MAT-file defining `ids` (list of test case IDs)|
| `m_imr_fd`       | Finite difference solver function              |
| `m_imr_spectral` | Spectral solver function                       |

## Outputs

- **Console output**:
  - L2 norm error for each finite-difference and spectral test.
  - Summary line: pass or fail with indices.
- **Exit code**:
  - `0` if all tests pass.
  - Error thrown (non-zero) if any test fails.

## Notes

- `num_tests` covers all combinations of radial, thermal, stress, and material parameters.
- Tests use reduced mesh sizes (`Nt=70, Mt=70` for FD; `Nt=12, Mt=12` for spectral) for speed.
