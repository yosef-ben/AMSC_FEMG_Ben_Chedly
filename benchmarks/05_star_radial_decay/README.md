# Benchmark 05: Four-Pointed Star, Radial Decay With Reaction

## Goal

Test a genuinely time-dependent parabolic problem on the four-pointed star,
including a linear reaction term `r u`. The benchmark has an exact solution, so
it can be checked both visually in ParaView and quantitatively with an L2 error.

## Domain

The domain is the four-pointed star in the plane:

```text
data/star_4.txt
```

The central vertex is located at `(0,0)`, and the four external vertices lie on
the coordinate axes.

## Problem

Find `u` on the metric graph such that

```text
u_t - Delta u + r u = 0
```

with homogeneous Dirichlet conditions at the four external vertices and the
natural Kirchhoff condition at the central vertex.

## Exact Solution

On each edge, let `s` be the local coordinate from the central vertex to the
external vertex. The exact solution is

```text
u(t,s) = exp(-(pi^2/4 + r)t) cos(pi s / 2).
```

Here `r = 1`. The function is zero at the external vertices and has zero
outgoing derivative at the center, so it satisfies the graph conditions.

## Numerical Parameters

```text
diffusion = 1
reaction  = 1
T         = 1
deltat    = 0.05
theta     = 0.5
n_cells   = 100 per edge
```

The temporal method is Crank-Nicolson.

## Reproduce

From the project root:

```bash
cmake --build build-release
./build-release/test_heat_star_radial_decay
./build-release/test_heat_star_radial_decay_time_convergence
python scripts/plot-time-convergence.py \
  output/convergence/star_radial_decay_time_convergence.csv \
  benchmarks/05_star_radial_decay/results/time_convergence.pdf
```

The run writes temporary visualization output to:

```text
output/visualization/star_radial_decay/
```

## Stored Results

The benchmark results selected for the report are stored in:

```text
benchmarks/05_star_radial_decay/results/
```

Open the following file in ParaView:

```text
benchmarks/05_star_radial_decay/results/solution.pvd
```

The temporal convergence files are:

```text
benchmarks/05_star_radial_decay/results/time_convergence.csv
benchmarks/05_star_radial_decay/results/time_convergence.pdf
```

## Result

The final L2 error against the exact radial decay solution is:

```text
3.86387e-04
```

The norm of the solution decreases monotonically from the initial profile to
the final time, as expected from diffusion plus the positive reaction term.

The temporal convergence test gives:

```text
dt     BE L2       BE rate  CN L2       CN rate
0.2    5.74135e-02 -        6.14347e-03 -
0.1    2.79361e-02 1.03926  1.53402e-03 2.00174
0.05   1.36640e-02 1.03175  3.83558e-04 1.99980
0.025  6.73801e-03 1.01998  9.61081e-05 1.99671
```

This confirms the expected first-order convergence of Backward Euler and
second-order convergence of Crank-Nicolson.

## Report Note

This is the first benchmark in which the solution genuinely evolves in time.
It also introduces the linear term `r u`, matching the structure of the
Hamiltonian/operator terms used in the elliptic graph formulation. The
temporal convergence table is the main quantitative validation of the
theta-method implementation.
