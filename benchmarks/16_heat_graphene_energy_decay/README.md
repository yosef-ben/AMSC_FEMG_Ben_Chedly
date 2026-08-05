# Benchmark 16: Energy Decay on a Graphene-Like Graph

## Goal

Study the dissipative behaviour of the heat equation on a nontrivial metric
graph. This benchmark is inspired by the time-dependent tests discussed by
Benzi for diffusion problems on quantum graphs.

Unlike the single-eigenmode tests, the initial datum is a combination of
several higher eigenmodes. Therefore the shape of the solution changes in
time: the highest-frequency components decay faster and the solution becomes
progressively smoother.

## Domain

The domain is the graphene-like graph:

```text
data/graphene_13.txt
```

It has 12 vertices, 13 unit edges, and 100 finite elements per edge.

## Problem

Solve

```text
u_t - Delta_Gamma u = 0
```

with homogeneous Neumann conditions at the external vertices and
Neumann-Kirchhoff conditions at the internal vertices.

The initial datum is built as a modal combination:

```text
u(0) = u_6 - 0.75 u_7 + 0.50 u_8 - 0.35 u_9,
```

where `u_k` denotes the k-th computed eigenmode of the graph Laplacian on the
same finite element space.

The associated eigenvalues are:

```text
lambda_6 = 3.58081209034
lambda_7 = 4.38664969359
lambda_8 = 4.38664969361
lambda_9 = 4.75036328508
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

## Diagnostics

At each time step we store:

```text
u_h(t)^T M u_h(t)
u_h(t)^T H u_h(t)
```

where `M` is the mass matrix and `H` is the stiffness matrix. For the
homogeneous heat equation both quantities are expected to decrease in time.

Since the initial condition is a linear combination of mass-normalized
eigenmodes, the modal reference is:

```text
||u(t)||_M^2 = sum_k c_k^2 exp(-2 lambda_k t),
E(t)         = sum_k lambda_k c_k^2 exp(-2 lambda_k t).
```

## Reproduce

From the project root:

```bash
cmake -S . -B build
cmake --build build --target test_heat_graphene_energy_decay
./build/test_heat_graphene_energy_decay
```

The temporary output is written to:

```text
output/visualization/graphene_energy_decay/
output/energy/graphene_energy_decay.csv
```

The stored benchmark output is in:

```text
benchmarks/16_heat_graphene_energy_decay/results/
```

Open in ParaView:

```text
benchmarks/16_heat_graphene_energy_decay/results/solution.pvd
```

## Result

The final values at `T = 1` are:

```text
u_h^T M u_h = 9.099538130993e-04
u_h^T H u_h = 3.370290165131e-03
```

Both the discrete mass norm and the discrete energy decay monotonically over
time. This confirms the dissipative character of the heat equation on the
graph and gives a more global parabolic diagnostic than a single snapshot.
