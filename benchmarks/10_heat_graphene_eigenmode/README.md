# Benchmark 10: Heat Equation on a Graphene Eigenmode

## Goal

Validate the parabolic solver on a nontrivial graph topology using an exact
time-dependent solution built from a Laplacian eigenfunction.

This benchmark is the parabolic extension of the spectral problem on the
graphene-like metric graph. If `phi_k` solves

```text
-Delta_Gamma phi_k = lambda_k phi_k,
```

then the heat equation has the exact solution

```text
u(t) = exp(-mu * lambda_k * t) * phi_k.
```

## Domain

The domain is the graphene-like graph:

```text
data/graphene_13.txt
```

It has 12 vertices, 13 unit edges, and 100 finite elements per edge.

## Problem

Solve

```text
u_t - mu Delta_Gamma u = 0
```

on the metric graph, with homogeneous Neumann conditions at the external
vertices and Kirchhoff conditions at the internal vertices.

The initial condition is the sixth eigenfunction of the graph Laplacian. This higher mode is still an exact modal solution, but its larger eigenvalue makes the temporal decay clearly visible in the snapshots:

```text
u(0) = phi_6.
```

The eigenvalue is

```text
lambda_6 = 3.580812090341592.
```

Therefore the exact solution is

```text
u(t) = exp(-3.580812090341592 * t) * phi_6,
```

since `mu = 1`.

## Numerical Parameters

```text
mu      = 1
T       = 1
deltat  = 0.01
theta   = 0.5
n_cells = 100 per edge
```

The theta method is Crank--Nicolson.

## Reproduce

From the project root:

```bash
cmake --build build --target test_heat_graphene_eigenmode
./build/test_heat_graphene_eigenmode
```

The temporary output is written to:

```text
output/visualization/graphene_eigenmode/
```

The stored benchmark output is in:

```text
benchmarks/10_heat_graphene_eigenmode/results/
```

The decay of the numerical L2 norm is stored in:

```text
benchmarks/10_heat_graphene_eigenmode/results/decay.csv
```

Open in ParaView:

```text
benchmarks/10_heat_graphene_eigenmode/results/solution.pvd
```

## Result

The final error at `T = 1` is

```text
Final L2 error = 1.0657e-05.
```

The solution keeps the shape of the eigenfunction and its amplitude decays
with the factor `exp(-lambda_6 t)`. This confirms the expected modal
behaviour of the heat equation on a metric graph.

The file `decay.csv` compares the numerical L2 norm with the exact decay
law. Since the eigenfunction is mass-normalized, the initial L2 norm is one.
