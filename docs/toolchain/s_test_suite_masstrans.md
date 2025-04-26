# `s_test_suite_masstrans.m` Documentation

## Overview
`s_test_suite_masstrans.m` validates the finite difference solver’s mass transfer implementation under various model settings.

## Major Sections

1. **Initialization**
   - Clears workspace.
   - Adds paths: `src/forward_solver/` and `tests/`.
   - Loads `file_ids.mat`.

2. **Preallocation**
   - Computes `num_tests`.
   - Initializes `errors_fd` and `failed_tests`.
   - Sets `threshold = 1e-5`, `count = 1`.

3. **Test Parameter Setup**
   - Defines `tvector` for simulation.
   - Fixes `masstrans = 1`, `vapor = 1`, `collapse = 1`.
   - Sets bubble geometry: `R0`, `Req`.
   - Defines model vectors: `radial`, `bubtherm`, `medtherm`, `stress`.

4. **Loop Over Test Cases**
   - Iterates over all combinations.
   - Builds `varin` parameters.
   - Loads reference `Rm` from test files.
   - Executes `m_imr_fd(...,'Nt',30,'Mt',70)`, computes L2 error.
   - Records failures above `threshold`.

5. **Summary**
   - Trims `failed_tests`.
   - Prints pass/fail summary (no automatic exit).

## Inputs

| Name            | Description                          |
|-----------------|--------------------------------------|
| `file_ids.mat`  | Test case identifiers                |
| `m_imr_fd`      | Finite difference solver function    |

## Outputs

- **Console output**: L2 norm errors per test and summary.
- **Error flags**: `failed_tests` indices where mass transfer failed.

## Notes

- Script does not exit on failure to allow further analysis.
- Uses moderate resolution (`Nt=30`, `Mt=70`) for faster execution.
