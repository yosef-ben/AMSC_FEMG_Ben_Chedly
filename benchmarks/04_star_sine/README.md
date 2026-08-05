# Benchmark 04: Four-Pointed Star, Sine Solution

## Goal

Validate the parabolic solver on a non-polynomial function over the
four-pointed star. This is the parabolic analogue of the sine-function
elliptic test in Abbate and Di Primio.

Unlike the constant and linear benchmarks, the exact solution is not contained
in the P1 finite element space. Therefore the numerical error is expected to be
nonzero and to decrease with mesh refinement.

## Domain

The domain is the four-pointed star in the plane:

```text
data/star_4.txt
```

For the spatial convergence test, the refined graph files are:

```text
data/star_4_10.txt
data/star_4_20.txt
data/star_4_40.txt
data/star_4_80.txt
```

## Problem

Find `u` on the metric graph such that

```text
u_t - Delta u = f
```

with homogeneous Dirichlet conditions at all original vertices of the star.

## Exact Solution

The prescribed stationary solution is

```text
u(t,x,y) = sin(2*pi*x).
```

On the two horizontal edges this gives a sine profile. On the two vertical
edges the solution is identically zero.

## Forcing Term

Since the exact solution is stationary,

```text
u_t = 0,
```

and the forcing term is chosen as

```text
f(x,y) = 4*pi^2*sin(2*pi*x).
```

## Numerical Parameters

For the visualization run:

```text
diffusion = 1
T         = 0.1
deltat    = 0.005
theta     = 1
n_cells   = 100 per edge
```

The temporal method is Backward Euler.

For the spatial convergence run:

```text
n_cells = 10, 20, 40, 80 per edge
deltat  = 0.025*h^2
T       = 0.01
theta   = 1
```

The small time step makes the spatial error dominant.

## Reproduce

From the project root:

```bash
cmake --build build
./build/test_heat_star_sine
./build/test_heat_star_sine_convergence
```

The visualization run writes temporary output to:

```text
output/visualization/star_sine/
```

The convergence run writes:

```text
output/convergence/star_sine_space_convergence.csv
```

## Stored Results

The benchmark results selected for the report are stored in:

```text
benchmarks/04_star_sine/results/
```

Open the following file in ParaView:

```text
benchmarks/04_star_sine/results/solution.pvd
```

The spatial convergence table is stored in:

```text
benchmarks/04_star_sine/results/space_convergence.csv
```

## Result

The visualization run gives the final L2 error:

```text
3.28946e-04
```

The spatial convergence run gives:

```text
h       L2 error     rate
0.1     3.25637e-02  -
0.05    8.20369e-03  1.98892
0.025   2.05486e-03  1.99724
0.0125  5.13960e-04  1.99931
```

This confirms the expected second-order L2 convergence for P1 finite elements.

## Report Note

This benchmark is the first one in which the exact solution is not reproduced
exactly by the finite element space. It is therefore useful for showing the
actual approximation behavior of the parabolic solver on a metric graph.
