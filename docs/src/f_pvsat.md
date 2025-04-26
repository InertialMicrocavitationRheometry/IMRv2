# `f_pvsat.m` Documentation

## Overview
`f_pvsat` returns the saturated water vapor pressure as a function of temperature using an empirical exponential fit.

## Major Sections
1. **Empirical Formula**: Applies `Pv = 1.17e11 * exp(-5200/T)`.

## Inputs

| Name | Description        | Units |
|------|--------------------|-------|
| `T`  | Temperature       | K     |

## Outputs

| Name | Description                   | Units |
|------|-------------------------------|-------|
| `Pv` | Saturated vapor pressure      | Pa    |

## Notes
- Source: A. Preston, 2004 (Thesis).
- Valid for typical biological temperature ranges (290–310 K).
