# Benchmark 09: Spectral Problem on a Binary Tree

## Goal

Compute the first eigenpairs of the graph Laplacian on a tree-shaped metric
graph. This benchmark introduces a branched topology with several internal
junctions and many leaves. It is useful both for comparison with spectral
results on quantum graphs and as a preparation for parabolic propagation on
tree-like networks.

As before, each eigenfunction can be used to manufacture an exact solution of
the heat equation:

```text
u(t) = exp(-mu * lambda_k * t) * phi_k.
```

## Domain

The graph is a binary tree with 16 vertices and 15 unit edges:

```text
data/tree_15.txt
```

Each edge is discretized with 100 elements, giving 1501 degrees of freedom.

## Problem

Find eigenpairs `(lambda, phi)` such that

```text
-Delta_Gamma phi = lambda phi
```

with homogeneous Neumann conditions at the leaves and Kirchhoff
continuity/flux balance at every internal junction.

The finite element discretization gives the generalized eigenvalue problem

```text
H phi = lambda M phi.
```

## Reproduce

From the project root:

```bash
cmake --build build --target test_spectral_tree
./build/test_spectral_tree data/tree_15.txt
```

The temporary output is written to:

```text
output/spectral/tree/
```

The stored benchmark output is in:

```text
benchmarks/09_spectral_tree/results/
```

## Visualization

Open in ParaView:

```text
benchmarks/09_spectral_tree/results/domain.vtp
benchmarks/09_spectral_tree/results/eigenmode_00.vtp
benchmarks/09_spectral_tree/results/eigenmode_01.vtp
...
```

The file `domain.vtp` contains the flat tree. The files `eigenmode_XX.vtp`
contain the eigenfunction lifted in the third coordinate. The report figures
use the flat metric graph in blue and the lifted eigenfunction in red. The
vertical scale is used only for visualization.

The generated summary figures are:

```text
benchmarks/09_spectral_tree/results/tree_eigenmodes_00_07.png
benchmarks/09_spectral_tree/results/tree_eigenmodes_08_15.png
```

## Result

The first eigenvalues are:

```text
lambda_0  = 1.51040215977e-12
lambda_1  = 0.115489236187
lambda_2  = 0.378816467632
lambda_3  = 0.378816467643
lambda_4  = 0.707400678706
lambda_5  = 2.46745183459
lambda_6  = 2.46745183459
lambda_7  = 2.46745183460
lambda_8  = 2.46745183460
lambda_9  = 2.46745183460
lambda_10 = 2.46745183461
lambda_11 = 5.29264401401
lambda_12 = 6.38158595362
lambda_13 = 6.38158595362
lambda_14 = 7.85034876253
lambda_15 = 9.87041617023
```

The zero eigenvalue corresponds to the constant eigenfunction. The repeated
eigenvalues are a consequence of the symmetries of the tree branches. Higher
eigenmodes become progressively more oscillatory along the graph.
