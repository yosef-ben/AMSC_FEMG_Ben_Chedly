# Benchmark 11: Heat Equation on a Tree Eigenmode

## Goal

Validate the parabolic solver on a branched tree topology using an exact
time-dependent solution built from a Laplacian eigenfunction.

If `phi_k` solves

```text
-Delta_Gamma phi_k = lambda_k phi_k,
```

then the heat equation has the exact solution

```text
u(t) = exp(-mu * lambda_k * t) * phi_k.
```

## Domain

The domain is the binary tree:

```text
data/tree_15.txt
```

It has 16 vertices, 15 unit edges, and 100 finite elements per edge.

## Problem

Solve

```text
u_t - mu Delta_Gamma u = 0
```

with homogeneous Neumann conditions at the leaves and Kirchhoff conditions at
the internal junctions. The initial condition is the fifth eigenfunction. This mode is chosen because its larger eigenvalue gives a visible temporal decay while preserving an exact solution:

```text
u(0) = phi_5.
```

The eigenvalue is

```text
lambda_5 = 2.467451834590349.
```

Therefore, with `mu = 1`, the exact solution is

```text
u(t) = exp(-2.467451834590349 * t) * phi_5.
```

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
cmake --build build --target test_heat_tree_eigenmode
./build/test_heat_tree_eigenmode
```

The stored benchmark output is in:

```text
benchmarks/11_heat_tree_eigenmode/results/
```

Open in ParaView:

```text
benchmarks/11_heat_tree_eigenmode/results/solution.pvd
```

The L2 decay comparison is stored in:

```text
benchmarks/11_heat_tree_eigenmode/results/decay.csv
```

## Result

The final error at `T = 1` is

```text
Final L2 error = 1.06164e-05.
```

The solution remains proportional to the selected eigenfunction and decays
according to the exact factor `exp(-lambda_5 t)`.
