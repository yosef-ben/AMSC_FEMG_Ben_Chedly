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

### Front speed: a diagnostic, not a verification

For `alpha=1` the measured front speeds exceed the Fisher-KPP asymptotic value
`2*sqrt(d*alpha)` by `1.9%`, `4.0%` and `5.2%` for increasing diffusion. The
figure `results/front_speeds.pdf` is kept here as a diagnostic, but it is
deliberately **not** used in the report, for three reasons established by
measurement.

The excess depends on the window over which the speed is fitted, by more than
a percentage point:

```text
d          whole run   t >= 10   t >= 14   t >= 16
1e-4            2.4%      1.9%      2.4%      2.7%
2e-4            5.0%      4.0%      4.5%      4.8%
4e-4            6.4%      5.2%      5.6%      5.8%
```

The three values quoted above are the smallest of each row, which is an
artefact of the window the diagnostic script happens to use.

The excess does not decrease when the fit is restricted to later times. A
finite-time correction to the asymptotic speed would approach it from below and
shrink; this one sits above and is flat, so it is not the transient one would
want to fit away.

The excess cannot be shown to converge, because the time step cannot be
refined. With the sharp single-node initial datum the consistent mass matrix
undershoots, and the undershoot **grows** under time refinement:

```text
dt        minimum concentration at d = 1e-4, alpha = 2
0.1       -1.5e-13
0.05      -5.4e-06
0.025     -9.1e-04
```

so the driver's physical-range guard fires at `dt = 0.05` and below. This is
the same loss of the discrete maximum principle that benchmark 23 establishes
on the connectome, and it is useful to have found it here as well: it is a
property of the consistent-mass P1 discretization, not of the graph. The
executable accepts the time step as a third argument so that this can be
reproduced.

What the report says instead is the qualitative statement the reference itself
makes: the fronts are symmetric and their speed grows with both `d` and
`alpha`, consistently with `2*sqrt(d*alpha)`.

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
results/time_step_study.pdf        (report figure: the profiles alone)
results/time_step_study_full.pdf   (both panels, kept as the record)
```

### What this study is, and what it is not

It is a **time-step sensitivity study**, not a convergence study. Weickenmeier
et al. describe the experiment in words only, in section 2.4 on page 268: they
report that `dt = 0.1` was regarded as sufficiently converged and that larger
steps produce a spurious increase in spreading. No data and no figure are
given, so their criterion cannot be reconstructed, and nothing here claims
their statement is wrong.

The diagnostic used here is the position `x_f` of the right front, the point
where the profile crosses `c = 0.5`. The finest run, `dt = 0.025`, is a
**numerical reference and not an exact solution**, so the difference
quantities below measure sensitivity to the step and not error. The formal temporal orders
are established elsewhere, against the exact logistic solution: `1.00` for
Backward Euler and `2.00` for the semi-implicit Crank-Nicolson scheme.

### How to read the figure

The report figure (`time_step_study.pdf`) draws all six runs and all 201 nodal
values of each; nothing is omitted. Colour encodes the time step, and the
cool-to-warm break falls exactly at the step above which no front survives.
The `dt = 0.4` profile is dashed because it coincides with `dt = 0.3` to
`4e-5`. On the plateau `|x| <= 0.5` all six runs agree to `3.3e-4`, so the
disagreement really is confined to the front.

The full variant (`time_step_study_full.pdf`, panel selection `--panel both`
of the plotting script, kept in the results and not in the report) adds a
second panel with three dimensionless differences from the finest run, all
evaluated on the profiles at the common final time `T = 19.2`, with
`Omega = (-1,1)` and `c_ref` the profile at `dt = 0.025`:

```text
e_inf(dt) = max_x |c_dt - c_ref|
e_2(dt)   = ( (1/|Omega|) * integral_Omega (c_dt - c_ref)^2 dx )^(1/2)
e_f(dt)   = |x_f(dt) - x_f(ref)|
```

`e_2` is a root mean square, that is the stored `l2_error` column divided by
`sqrt(|Omega|) = sqrt(2)`, so that it is dimensionless and comparable with
`e_inf`; `e_f` is a displacement in a domain of half-width one, hence also
dimensionless. The values drawn are

```text
dt        e_inf      e_2        e_f
0.025     0          0          0        (reference, absent from the log axis)
0.05      0.434688   0.106268   0.037930
0.1       0.898898   0.288002   0.125671
0.2       0.998928   0.534894   does not exist
0.3       0.999992   0.542931   does not exist
0.4       1.000000   0.542938   does not exist
```

`e_f` is defined only where the profile still crosses `c = 0.5`, so its curve
stops at `dt = 0.1`; the shaded region of that panel records that the quantity
does not exist beyond it, and no censored marker is drawn. The two norms of `c` remain well
defined for every step and saturate against the uniformly filled state.

No slopes are drawn on the figure. They were removed deliberately: on a figure
whose reference is another numerical run, an annotated slope invites the reader
to take it for a convergence order. The numbers are recorded below instead.

The figure follows the visual conventions of the reference's line figures: full
frame, no grid, sparse inward ticks, bold labels with the axis name after the
last tick, direct coloured series labels and no legend box. It is drawn at
`7.2` by `3.05` inches so that it stays legible when set at text width.

### The two-point slopes, recorded here and not shown

```text
measure            dt = 0.05 -> 0.1     dt = 0.1 -> 0.2
L2 error                     1.438                0.89
L-infinity error             1.048                0.15
front position               1.728                 --
```

Only two of the six steps still produce a solution with a front, so these are
two-point slopes and not fitted orders. They are kept in this record because
they are part of what the experiment produced, but they are not drawn on the
figure and should not be quoted as a recovery of the first order of Backward
Euler, for three reasons: the second column is
saturation rather than a rate; at `dt = 0.1` the L-infinity error already sits
at `0.899` of its own ceiling of `1`; and the errors are measured against the
finest run rather than an exact solution, which biases the estimate. A genuine
order would need at least two more steps below `0.1` with a finer reference,
and a longer domain so that the front never reaches `x = 1` within `T = 19.2`.

### Conventions in the stored CSV

`l2_error` is `sqrt(int_{-1}^{1} (c - c_ref)^2)` by the trapezoid rule on the
201 nodes, and `max_error` is `max_x |c - c_ref|`, both against `dt = 0.025`.
`mean_concentration` is the trapezoidal average over `[-1,1]` divided by two,
not the plain sample mean; the sample mean gives `0.6822` where the stored
value is `0.6856`. It is not plotted, because it equals the front position to
within `0.0044` at every tested step. `front_position = 1` for
`dt = 0.2, 0.3, 0.4` is a sentinel written when no `c = 0.5` crossing exists:
the true minima there are `0.93200`, `0.99996` and `0.99999994`, so the front
has left the interval entirely and those three points are lower bounds rather
than measured positions.

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
