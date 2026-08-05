# Benchmark 02: Four-Pointed Star, Constant Solution

## Goal

Verify that the parabolic solver preserves constant functions on a nontrivial
metric graph, following the first numerical experiment of Abbate and Di Primio
but in the parabolic setting.

## Domain

The domain is a four-pointed star in the plane with 5 original vertices and 4
unit-length edges. The central vertex is connected to four external vertices
lying on the coordinate axes.

Input graph:

```text
data/star_4.txt
```

## Problem

Find `u` on the metric graph such that

```text
u_t - Delta u = 0
```

with natural Kirchhoff conditions at the graph vertices.

## Initial Condition

```text
u_0 = 5
```

## Exact Solution

```text
u(t, x, y) = 5
```

for all times and on every edge.

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
./build/test_heat_star_constant
```

The run writes temporary visualization output to:

```text
output/visualization/star_constant/
```

## Stored Results

The benchmark results selected for the report are stored in:

```text
benchmarks/02_star_constant/results/
```

Open the following file in ParaView:

```text
benchmarks/02_star_constant/results/solution.pvd
```

## Result

The final L2 error against the exact constant solution is:

```text
4.25689e-12
```

The error is at roundoff level, confirming that the method preserves constant
solutions on the four-pointed star.

## Report Note

This benchmark is the parabolic analogue of the constant-function elliptic test
in Abbate and Di Primio. It checks continuity at the central vertex and the
correct assembly of mass and stiffness contributions on a nontrivial graph.
