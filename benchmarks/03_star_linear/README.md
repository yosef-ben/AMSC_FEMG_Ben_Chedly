# Benchmark 03: Four-Pointed Star, Linear Solution

## Goal

Verify that the parabolic solver preserves a linear function on the
four-pointed star when compatible Dirichlet values are imposed at the external
vertices. This is the parabolic analogue of the linear-function elliptic test
in Abbate and Di Primio.

## Domain

The domain is the four-pointed star in the plane:

```text
data/star_4.txt
```

The central vertex is located at `(0,0)`, and the four external vertices lie on
the coordinate axes.

## Problem

Find `u` on the metric graph such that

```text
u_t - Delta u = 0
```

with Kirchhoff condition at the central vertex and Dirichlet conditions at the
four external vertices.

## Initial Condition

```text
u_0(x,y) = x + y
```

## Boundary Values

The external vertices have the following values:

```text
u( 1,  0) =  1
u( 0,  1) =  1
u(-1,  0) = -1
u( 0, -1) = -1
```

## Exact Solution

```text
u(t,x,y) = x + y
```

for all times. The function is linear on each edge and satisfies Kirchhoff
continuity at the central vertex.

## Numerical Parameters

```text
diffusion = 1
T         = 1
deltat    = 0.05
theta     = 1
n_cells   = 100 per edge
```

The temporal method is Backward Euler.

## Reproduce

From the project root:

```bash
cmake --build build
./build/test_heat_star_linear
```

The run writes temporary visualization output to:

```text
output/visualization/star_linear/
```

## Stored Results

The benchmark results selected for the report are stored in:

```text
benchmarks/03_star_linear/results/
```

Open the following file in ParaView:

```text
benchmarks/03_star_linear/results/solution.pvd
```

## Result

The final L2 error against the exact linear solution is:

```text
4.348e-14
```

The error is at roundoff level, confirming that the method preserves compatible
linear solutions on the four-pointed star.

## Report Note

This benchmark validates the newly added Dirichlet support and checks that the
discrete solution remains globally continuous at the central vertex.
