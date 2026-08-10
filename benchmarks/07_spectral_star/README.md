# Benchmark 07: Spectral Problem on the Four-Pointed Star

## Goal

Compute the first eigenpairs of the graph Laplacian on the four-pointed star.
This is the first step of the spectral pipeline: once an eigenfunction
`phi_k` is known, the heat equation has the exact solution

```text
u(t) = exp(-mu * lambda_k * t) * phi_k.
```

Therefore this benchmark will be used as the stationary spectral reference for
the following parabolic tests.

## Domain

The graph is the four-pointed star:

```text
data/star_4.txt
```

It has one central vertex and four external vertices. Each edge has length 1
and is discretized with 100 elements.

## Problem

Find eigenpairs `(lambda, phi)` such that

```text
-Delta_Gamma phi = lambda phi
```

with homogeneous Neumann conditions at the four external vertices and
Kirchhoff continuity/flux balance at the central vertex.

In matrix form, the finite element discretization gives the generalized
eigenvalue problem

```text
H phi = lambda M phi,
```

where `H` is the stiffness matrix and `M` is the mass matrix.

## Reproduce

From the project root:

```bash
cmake --build build-release --target test_spectral_star
./build-release/test_spectral_star data/star_4.txt
```

The temporary output is written to:

```text
output/spectral/star/
```

The stored benchmark output is in:

```text
benchmarks/07_spectral_star/results/
```

## Visualization

Open in ParaView:

```text
benchmarks/07_spectral_star/results/domain.vtp
benchmarks/07_spectral_star/results/eigenmode_00.vtp
benchmarks/07_spectral_star/results/eigenmode_01.vtp
...
```

The file `domain.vtp` contains the flat graph. The files `eigenmode_XX.vtp`
contain the eigenfunction lifted in the third coordinate, so that the
oscillations are visible without using the Tube filter. The scalar field is
named `phi`.

## Result

The first eigenvalues are:

```text
lambda_0 = 4.87053651966e-12
lambda_1 = 2.46745183456
lambda_2 = 2.46745183460
lambda_3 = 2.46745183460
lambda_4 = 9.87041617021
lambda_5 = 22.2107196526
lambda_6 = 22.2107196526
lambda_7 = 22.2107196526
lambda_8 = 39.4914071916
lambda_9 = 61.7167427111
```

The zero eigenvalue corresponds to the constant mode. The first nonzero
eigenvalue is close to `(pi/2)^2` and has multiplicity three, as expected from
the symmetry of the four-pointed star. The next simple eigenvalue is close to
`pi^2`.
