# Benchmark 13: Combinatorial and Metric Laplacian Spectra

## Goal

Compare the spectrum of the combinatorial graph Laplacian with the spectrum
of the finite element Laplacian on the corresponding metric graph.

This benchmark is stationary, but it completes the spectral validation before
using the metric eigenpairs to construct exact parabolic solutions.

## Problems Compared

For the original graph, the combinatorial Laplacian is

```text
L = D - A,
```

where `D` is the degree matrix and `A` is the adjacency matrix.

For the metric graph, the finite element discretization gives the generalized
eigenvalue problem

```text
H phi = lambda M phi.
```

The comparison is meaningful between the `N` combinatorial eigenvalues and
the first `N` metric eigenvalues, where `N` is the number of original vertices.

## Graphs

The benchmark is run on:

```text
data/star_4.txt
data/graphene_13.txt
data/tree_15.txt
```

## Reproduce

From the project root:

```bash
cmake --build build-release --target test_spectral_comparison
./build-release/test_spectral_comparison
```

The stored outputs are:

```text
benchmarks/13_spectral_comparison/results/star/spectral_comparison.csv
benchmarks/13_spectral_comparison/results/graphene/spectral_comparison.csv
benchmarks/13_spectral_comparison/results/tree/spectral_comparison.csv
```

## Result

For the star graph, the comparison is:

```text
index,combinatorial_laplacian,metric_laplacian
0,0,4.870536519657767e-12
1,1,2.467451834559538
2,1,2.467451834595166
3,1,2.467451834597808
4,5,9.870416170214265
```

For the graphene-like graph, the first values are:

```text
index,combinatorial_laplacian,metric_laplacian
0,0,5.850707553677823e-12
1,0.1577643206756957,0.1570615633371609
2,1,0.925651336854093
3,1,1.096632732782526
4,1,1.096632732783475
5,1.493058083849815,1.560814525324424
```

For the tree, the first values are:

```text
index,combinatorial_laplacian,metric_laplacian
0,0,1.510402159771394e-12
1,0.09678807408844715,0.1154892361869879
2,0.2679491924311225,0.3788164676323109
3,0.2679491924311227,0.3788164676427527
4,0.4964955010604721,0.7074006787056091
```

The comparison confirms that the metric Laplacian preserves qualitative
features of the combinatorial spectrum, such as zero eigenvalues and
multiplicities induced by graph symmetries, while the numerical values differ
because the metric operator also depends on edge lengths and finite element
discretization.
