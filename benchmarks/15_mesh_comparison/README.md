# Benchmark 15: Mesh Comparison on a Nonuniform Tree

## Goal

Compare different discretization strategies on a metric tree with nonuniform
edge lengths. This benchmark adapts the mesh-comparison study to the current
code base and prepares the parabolic extension.

The comparison is made on the spectrum of the metric graph Laplacian. This is
relevant for the heat equation because each eigenmode evolves as

```text
u_k(t) = exp(-lambda_k t) phi_k.
```

Therefore, errors in the computed eigenvalues directly affect the modal
decay rates of the parabolic problem.

## Domain

The topology is the same binary tree used in the previous tests, but the edge
lengths decrease with the bifurcation level:

```text
level 0: length 1
level 1: length 1
level 2: length 0.5
level 3: length 0.25
```

## Meshes

Three meshes are compared:

```text
N-type  : same number of cells on every edge
h-type  : approximately uniform spatial step, comparable number of DoFs
h-min   : approximately uniform spatial step equal to the shortest N-type step
```

The generated graph files are stored in:

```text
benchmarks/15_mesh_comparison/results/tree_nonuniform_n_type.txt
benchmarks/15_mesh_comparison/results/tree_nonuniform_h_type.txt
benchmarks/15_mesh_comparison/results/tree_nonuniform_h_min.txt
```

The degrees of freedom are:

```text
mesh    dofs
N-type  901
h-type  900
h-min   1681
```

## Reproduce

From the project root:

```bash
cmake --build build-release --target test_spectral_mesh_comparison
./build-release/test_spectral_mesh_comparison
```

The stored outputs are:

```text
benchmarks/15_mesh_comparison/results/tree_mesh_summary.csv
benchmarks/15_mesh_comparison/results/tree_mesh_spectra.csv
```

## Result

The first eigenvalues are:

```text
index,N-type,h-type,h-min
0,0,0,0
1,0.3831895554,0.3831894461,0.3831889131
2,1.1826457938,1.1826285467,1.1826235345
3,2.8296061495,2.8296267581,2.8295976180
4,2.8296061495,2.8296267581,2.8295976180
5,4.4173588897,4.4170604460,4.4169907649
6,7.1989703713,7.1983230569,7.1981374385
```

For the lowest part of the spectrum, the N-type and h-type meshes agree well
when the number of degrees of freedom is comparable. The h-min mesh provides a
finer reference. The largest eigenvalues are expected to be more sensitive to
the mesh distribution, because their eigenfunctions are more oscillatory.

## Parabolic Interpretation

In the parabolic problem, these eigenvalues determine the decay factors

```text
exp(-lambda_k t).
```

Thus this benchmark justifies using h-type meshes on nonuniform graphs before
running time-dependent simulations on tree-like or connectome-like networks.
