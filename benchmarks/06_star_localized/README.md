# Benchmark 06: Eigenmode Relaxation on the Four-Pointed Star

## Goal

Validate the parabolic solver on a time-dependent exact solution built from an
eigenmode of the graph Laplacian. The initial condition is asymmetric: one
branch is hotter than the other three. The nonconstant part then decays in
time and the solution relaxes toward a constant state.

## Domain

The domain is the four-pointed star in the plane:

```text
data/star_4.txt
```

## Problem

Find `u` on the metric graph such that

```text
u_t - mu Delta u = 0
```

with the Kirchhoff condition at the central vertex and natural Neumann
conditions at the four external vertices.

## Exact Solution

On each edge, let `s` be the local coordinate from the central vertex to the
external vertex. The exact solution is

```text
u_i(t,s) = C + A a_i exp(-mu*pi^2*t/4) sin(pi*s/2),
```

with

```text
C = 4
A = 2
(a_1, a_2, a_3, a_4) = (3, -1, -1, -1).
```

The condition `a_1 + a_2 + a_3 + a_4 = 0` gives the Kirchhoff condition at the
central vertex.

## Numerical Parameters

```text
diffusion = 5
reaction  = 0
T         = 1
deltat    = 0.05
theta     = 0.5
n_cells   = 100 per edge
```

The temporal method is Crank-Nicolson.

## Reproduce

From the project root:

```bash
cmake --build build
./build/test_heat_star_localized
```

The run writes temporary visualization output to:

```text
output/visualization/star_localized/
```

## Stored Results

The benchmark results selected for the report are stored in:

```text
benchmarks/06_star_localized/results/
```

Open the following file in ParaView:

```text
benchmarks/06_star_localized/results/solution.pvd
```

## Result

The final L2 error against the exact eigenmode relaxation is:

```text
7.30537e-06
```

The solution relaxes toward the constant value `C = 4`.
