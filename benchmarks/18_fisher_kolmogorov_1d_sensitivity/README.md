# Benchmark 18: One-Dimensional Fisher-Kolmogorov Sensitivity

## Goal

Validate the complete nonlinear diffusion-reaction solver before applying it
to a brain connectome. The experiment follows the one-dimensional sensitivity
test of Weickenmeier et al. and exercises diffusion, nonlinear reaction,
finite-element assembly, and nonlinear time integration simultaneously.

## Problem

On the interval `[-1,1]`, solve

```text
c_t - d c_xx = alpha c (1-c)
```

with homogeneous Neumann boundary conditions. The initial concentration is
zero at every mesh node except the central node, where `c(0,0)=0.1`.

## Numerical Parameters

The interval is discretized with 200 P1 finite elements. The final time is
`T=20`, the time step is `dt=0.1`, and the nonlinear system at each step is
solved by Newton's method applied to Backward Euler.

The nine cases combine

```text
d     = 1e-4, 2e-4, 4e-4
alpha = 1, 2, 4
```

and the solution is sampled throughout the complete time interval.

## Reproduce

Run the commands in `commands.txt` from the project root.

## Results

The generated profiles and report-ready figure are stored in `results/`.
The expected behavior is a symmetric travelling front: increasing `alpha`
accelerates local conversion, while increasing `d` broadens and accelerates
spatial propagation. Concentrations should remain in the physical interval
`[0,1]`.

## Verification

The uniform logistic test verifies both nonlinear time integrators against the
exact scalar solution. The measured rates are approximately `2.00` for the
Corti semi-implicit scheme and `1.00` for Backward Euler.

Across the full sensitivity dataset, the maximum left-right symmetry error is
`3.33e-16`. The concentration range is
`[-6.34e-10, 1.00]`; the tiny negative value is roundoff-level numerical
noise. The final mean concentration increases consistently with either
diffusion or reaction strength, as expected.

The report-ready sensitivity figure and front-speed diagnostics are:

```text
results/sensitivity.pdf
results/front_speeds.csv
```

For `alpha=1`, the measured front speeds differ from the Fisher--KPP
asymptotic value `2*sqrt(d*alpha)` by approximately `1.9%`, `4.0%`, and
`5.2%` for increasing diffusion. This provides an independent check of the
coupled diffusion-reaction dynamics.

## Nonlinear time-step sensitivity

A second experiment follows the preliminary time-step study reported by
Weickenmeier et al. It fixes

```text
d = 2e-4
alpha = 2
dt = 0.025, 0.05, 0.1, 0.2, 0.3, 0.4
```

and compares all final profiles at the common time `T=19.2`. This time is
used instead of `T=20` because it is exactly divisible by every tested time
step, including `dt=0.3`. The `dt=0.025` solution is used only as the
finest numerical reference.

The experiment confirms that this nonlinear problem is highly sensitive to
the time step. The measured right-front positions are `0.6889`, `0.7268`,
and `0.8145` for `dt=0.025`, `0.05`, and `0.1`. For `dt>=0.2`, the
computed solution has already spread across almost the complete interval.
Thus large time steps preserve boundedness but introduce severe artificial
spreading, consistently with the warning in the reference study.

The corresponding files are:

```text
results/time_step_profiles.csv
results/time_step_study.csv
results/time_step_study.pdf
```

## Report figure

The numerical solution is sampled every `0.5` time units for visualization.
The report figure displays the complete space-time field `c(x,t)`: columns
correspond to increasing `alpha`, rows to increasing `d`, and every panel
uses the same color range `[0,1]`.

Suggested report caption:

> Space-time evolution of the misfolded-protein concentration for the
> one-dimensional Fisher-Kolmogorov problem. The growth coefficient increases
> from left to right, while the diffusion coefficient increases from top to
> bottom. Increasing alpha accelerates local conversion and indirectly the
> propagation front; increasing d directly enhances spatial spreading. All
> panels show the numerical FEM solution and use the same concentration scale.

## Literature verification

Checked against Weickenmeier, Jucker, Goriely and Kuhl, *A physics-based model
explains the prion-like features of neurodegeneration*, J. Mech. Phys. Solids
124 (2019) 264-281, sections 2.3 and 2.4.

Every element of the setup is taken from the paper and matches:

| Quantity | Paper (section 2.4) | This benchmark |
|---|---|---|
| domain | `B = {x : -1 <= x <= 1}` | same |
| elements | `n_el = 200` linear | 200 P1 |
| final time and steps | `T = 20`, `n_step = 200`, `dt = 0.1` | same |
| growth | `[1a, 2a, 4a]` with `a = 1` | `alpha = 1, 2, 4` |
| spreading | `[1d, 2d, 4d]` with `d = 0.0001` | `d = 1e-4, 2e-4, 4e-4` |
| initial datum | `c_0 = 0` except `c_0 = 0.1` at the centre node | same |
| boundary condition | homogeneous Neumann | same |
| time integration | implicit Euler backward | Backward Euler |
| nonlinear solver | Newton-Raphson on the consistently linearised residual, equations (10) and (11) | same |
| time-step study | `[dt/4, dt/2, dt, 2dt, 3dt, 4dt]` at `dt = 0.1`, `alpha = 2`, `d = 0.0002` | `dt = 0.025 ... 0.4`, same `alpha` and `d` |

The report figure `results/sensitivity.pdf` uses the same rendering as figure 3
of the paper: filled contours of the activation time in the `(x, c)` plane,
growth increasing across columns and spreading increasing down rows.

The paper reports the time-step study only qualitatively, as "a spurious
increase in spreading for larger time step sizes". The measured front positions
`0.6889`, `0.7268` and `0.8145` for `dt = 0.025, 0.05, 0.1`, against full
spreading for `dt >= 0.2`, make that statement quantitative.

`results/front_speeds.pdf` is an addition, not a reproduction: the paper gives
no front-speed check. Comparing the measured front against the Fisher-KPP
asymptotic speed `2 sqrt(d alpha)` is an independent verification of the
coupled dynamics, and the left panel confirms that the front becomes linear in
time before the domain boundary is reached.

All numbers in this file were regenerated from source in this working copy and
reproduce exactly.
