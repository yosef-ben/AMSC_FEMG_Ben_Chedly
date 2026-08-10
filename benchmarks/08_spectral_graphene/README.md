# Benchmark 08: Spectral Problem on a Graphene-Like Graph

## Goal

Compute the first eigenpairs of the graph Laplacian on a graphene-like metric
graph. This benchmark introduces a more complex topology than the
four-pointed star while keeping the geometry easy to inspect in ParaView.

The spectral problem is the stationary reference for the parabolic extension:
if `phi_k` is an eigenfunction associated with `lambda_k`, then the heat
equation admits the exact solution

```text
u(t) = exp(-mu * lambda_k * t) * phi_k.
```

## Domain

The graph is a planar graphene-like network with 12 vertices and 13 unit
edges:

```text
data/graphene_13.txt
```

Each edge is discretized with 100 elements, giving 1299 degrees of freedom.

## Problem

Find eigenpairs `(lambda, phi)` such that

```text
-Delta_Gamma phi = lambda phi
```

with homogeneous Neumann conditions at the external vertices and
Kirchhoff continuity/flux balance at the internal junctions.

The finite element matrices satisfy the generalized eigenvalue problem

```text
H phi = lambda M phi.
```

## Reproduce

From the project root:

```bash
cmake --build build-release --target test_spectral_graphene
./build-release/test_spectral_graphene data/graphene_13.txt
```

The temporary output is written to:

```text
output/spectral/graphene/
```

The stored benchmark output is in:

```text
benchmarks/08_spectral_graphene/results/
```

## Visualization

Open in ParaView:

```text
benchmarks/08_spectral_graphene/results/domain.vtp
benchmarks/08_spectral_graphene/results/eigenmode_00.vtp
benchmarks/08_spectral_graphene/results/eigenmode_01.vtp
...
```

Use `domain.vtp` as the flat red domain and color each `eigenmode_XX.vtp` by
the scalar field `phi`. The eigenfunction is lifted in the third coordinate,
so the oscillatory shape is visible directly.

## Result

The first eigenvalues are:

```text
lambda_0 = 5.85070755368e-12
lambda_1 = 0.157061563337
lambda_2 = 0.925651336854
lambda_3 = 1.09663273278
lambda_4 = 1.09663273278
lambda_5 = 1.56081452532
lambda_6 = 3.58081209034
lambda_7 = 4.38665119145
lambda_8 = 4.38665119146
lambda_9 = 4.75035862488
```

The zero eigenvalue corresponds to the constant eigenfunction. The repeated
eigenvalues reflect the symmetries of the graph topology.
