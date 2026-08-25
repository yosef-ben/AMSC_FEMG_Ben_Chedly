# Solver verification

Repository-internal correctness checks for the Fisher-Kolmogorov metric-graph
solver of `FEMG`. They answer one question, independently of the connectome
data and of any biological comparison: is the solver solving the intended
equations correctly? Nothing in this directory feeds the report.

The semi-discrete system the library implements is

```text
M c' + H c = R(c),    R_i(c) = sum_e int_e alpha_h c_h (1 - c_h) phi_i ds
```

with P1 elements on every edge, Kirchhoff coupling at the vertices, the
consistent mass (3-point Gauss on the quartic reaction integrand, exact) or
the row-sum lumped mass with the vertex rule, and two time schemes: backward
Euler with a damped, safeguarded Newton iteration, and the semi-implicit
Crank-Nicolson of Corti et al., in which the factor `(1 - c)` is
extrapolated to second order.

## Run

```bash
cmake --build build-release --target verify_fisher_kolmogorov
./build-release/verify_fisher_kolmogorov verification/results
```

One PASS/FAIL line per check, a machine-readable copy in
`results/verification_summary.csv`, exit code 0 only if every check passes.
The small graph files the tests use are written into the results directory
at run time. The whole suite runs in about half a minute.

## The tests

**1. Exact solution (interval and star).** With `alpha = 0` and unit edges,
`c(xi, t) = 1/2 + (1/4) e^{-pi^2 t} cos(pi xi)` solves the metric-graph heat
equation exactly: it has zero flux at both ends of each edge, so the
Kirchhoff sum at the centre of the symmetric three-edge star vanishes and no
forcing term is needed (the library has no source hook, so the test is built
on exact solutions rather than manufactured forcings). Spatial refinement at
`h = 1/4 .. 1/32` under a tiny Crank-Nicolson step measures the L2 rate
(expected 2) and the H1-seminorm rate (expected 1) by 3-point Gauss
quadrature against the exact mode; the absolute L2 error is checked against
the sharp interpolation bound `(h^2/pi^2) |u|_{H2} = 1.42e-4` at `h = 1/32`
(times `sqrt(3)` on the star), with 25 percent headroom because the finite
element solution is not the interpolant. Temporal refinement at fixed
`h = 1/64` measures rate 1 for backward Euler against the exact solution and
rate 2 for Crank-Nicolson against a small-step reference of itself, which
removes the spatial floor.

**2. Pure diffusion against the eigendecomposition.** On a graph with a
cycle, unequal lengths and diffusivities and a degree-3 vertex, the exact
solution of the semi-discrete system is
`c(t) = sum_k v_k e^{-lambda_k t} (v_k' M c_0)` with `H v_k = lambda_k M
v_k`, computed independently with a dense generalized eigensolver. Both
schemes and both masses are compared against it under step refinement:
rates 1 (backward Euler) and 2 (Crank-Nicolson), absolute M-norm errors at
the smallest step of order `1e-3` and `3e-6`.

**3. Pure reaction.** With the diffusivities set to zero the behaviour
depends on the mass:

- lumped mass with the vertex rule: the system decouples into one logistic
  equation per degree of freedom, with the exact solution
  `c(t) = c0 e^{at} / (1 + c0 (e^{at} - 1))`; backward Euler converges to it
  at rate 1;
- consistent mass, uniform state: the reaction vector is `alpha c (1-c)`
  times the mass row sums, so the semi-discrete solution is the scalar
  logistic exactly, and backward Euler again converges at rate 1;
- consistent mass, non-uniform state: the reaction couples neighbouring
  nodes through `M^{-1}`, so the semi-discrete flow is *not* nodewise
  logistic (the test measures the gap, 0.106 at the final time, so the
  distinction is real); the reference is a 40000-step Runge-Kutta
  integration of `M c' = R(c)` built from the same assembled reaction
  vector, and backward Euler converges to it at rate 1.

**4. Newton Jacobian.** For the backward Euler residual
`F(c) = M (c - c_old) + dt H c - dt R(c)` the assembled Jacobian
`J = M + dt H - dt W(c)` with `W_ij = int alpha_h (1 - 2 c_h) phi_i phi_j`
is checked against central finite differences of `F` in five random
admissible states and directions, over increments `1e-3 .. 1e-9`. The
directional relative error reaches round-off (about `3e-13`); the bound is
`1e-9`. Both masses.

**5. Space/time refinement consistency.** On the asymmetric graph with
nonuniform `alpha` and diffusivities, vertex values at `T = 8`: refining `h`
at `dt = 8/3200` gives Richardson differences shrinking at the P1 rate
(measured 1.77 over three pairs, slack 0.4, the coarsest level being one
cell on the shortest edge); refining `dt` at the finest mesh gives the
backward Euler rate 1.

**6. Invariants.**

- `alpha = 0`: `1' M c` is conserved by both schemes and both masses
  (row sums of `H` vanish); measured drift at round-off, bound `1e-12`
  relative;
- a uniform state with `alpha = 0` is stationary, and the reaction
  equilibria `c = 0` and `c = 1` are stationary with `alpha > 0`
  (`R(0) = R(1) = 0` exactly); measured deviation `0`;
- on the symmetric star with symmetric data the three edges carry the same
  solution; measured asymmetry `0`.

## What is deliberately not claimed

The test suite verifies the discretization the library implements, not the
positivity of the consistent-mass scheme (the loss of positivity in the
reaction-dominated regime is a documented property, benchmark 23) and not
the modelling choices of the applications. Test 1 needs no forcing because
the initial data are chosen among exact solutions; a manufactured forcing
would require a source-term hook the library does not have, and adding one
for testing alone was judged a larger intervention than the tests justify.

## Regression use

`verify_fisher_kolmogorov` fails (nonzero exit) if the mass or stiffness
assembly, the reaction assembly, the Newton Jacobian, either time scheme or
the vertex coupling changes incorrectly. Run it after touching
`FEMG/src/fisher_kolmogorov_problem.cpp` or the graph infrastructure,
together with the byte-identity regressions of the stored benchmark runs
(`benchmarks/*/commands.txt`) and `scripts/verify-figures.py`.

The one library change made for this suite is the exposure of
`assemble_reaction_weight_matrix` in the public interface of
`fisher_kolmogorov_problem` (declaration moved, no behaviour change); after
the change the stored lumped fornari83 staging run and the semi-implicit
corti83 run were re-run and reproduce their stored results byte for byte.
