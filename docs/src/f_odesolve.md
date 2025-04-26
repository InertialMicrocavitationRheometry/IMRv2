# `f_odesolve.m` Documentation

## Overview
`f_odesolve` is a wrapper that selects and executes MATLAB's ODE solvers (`ode15s`, `ode23tb`, `ode45`) with user-defined time stepping and tolerance options.

## Major Sections
1. **Options Configuration**: Sets `RelTol`, `AbsTol`, and `MaxStep` based on `divisions`.
2. **Solver Selection**: Routes to the correct solver based on `method` code:
   - 15 → `ode15s`
   - 23 → `ode23tb`
   - 45 → `ode45`
3. **Execution**: Calls the selected solver and returns results.

## Inputs

| Name       | Description                                  | Units   |
|------------|----------------------------------------------|---------|
| `bubble`   | Function handle for ODE system               | –       |
| `init`     | Initial condition vector                     | –       |
| `method`   | Solver choice (15, 23, or 45)                | –       |
| `divisions`| Maximum number of time steps (0 = auto)      | –       |
| `tspan`    | Two-element time span `[t0 tf]`              | –       |
| `tfin`     | Final time (for step size calculation)       | –       |

## Outputs

| Name | Description              | Units   |
|------|--------------------------|---------|
| `t`  | Time vector              | –       |
| `X`  | Solution matrix (states) | –       |

## Notes
- Throws an error for unsupported `method` codes.
- `MaxStep` is set to `tfin/divisions` when `divisions > 0`.
