# Benchmark 14: Parabolic Solver Complexity on the Graphene Graph

## Goal

Measure the computational cost of the parabolic solver as the number of
degrees of freedom increases. This benchmark adapts the computational
complexity study to the time-dependent setting.

The measured quantities are:

```text
init_seconds
assembly_seconds
factorization_seconds
time_stepping_seconds
total_solve_seconds
total_seconds
```

The output is disabled during this benchmark, so the timings refer to the
numerical solver and not to VTK file generation.

## Domain

The graph is the graphene-like topology. The same graph is discretized with
an increasing number of elements per edge:

```text
n_cells = 100, 200, 400, 800, 1600
```

The corresponding number of degrees of freedom is:

```text
1299, 2599, 5199, 10399, 20799
```

## Problem

Solve the homogeneous heat equation

```text
u_t - Delta_Gamma u = 0
```

using the default smooth initial condition. The final time and time step are:

```text
T      = 0.1
deltat = 0.01
```

Therefore each run performs 10 time steps. The theta method is
Crank--Nicolson.

## Reproduce

From the project root:

```bash
cmake --build build-release --target test_heat_graphene_complexity
./build-release/test_heat_graphene_complexity
```

The stored output is:

```text
benchmarks/14_parabolic_complexity/results/graphene_parabolic_complexity.csv
```

Each case is repeated three times and the reported values are averages.

## Result

```text
n_cells  n_dofs  total_seconds
100      1299    0.00343220
200      2599    0.00654443
400      5199    0.01281429
800      10399   0.02822035
1600     20799   0.05600564
```

The measured cost grows approximately linearly with the number of degrees of
freedom in this range. This is consistent with the sparse one-dimensional
structure of the metric graph problem and with the fact that the number of
time steps is fixed.
