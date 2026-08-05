# Benchmark 12: Time Convergence on a Graphene Eigenmode

## Goal

Verify the temporal order of the theta method on a nontrivial graph topology
using an exact modal solution of the heat equation.

The exact solution is

```text
u(t) = exp(-lambda_1 * t) * phi_1,
```

where `phi_1` is the first nonconstant eigenfunction of the graphene-like
graph and

```text
lambda_1 = 0.1570615633371609.
```

## Domain

The domain is the graphene-like graph:

```text
data/graphene_13.txt
```

with 100 finite elements per edge.

## Problem

Solve

```text
u_t - Delta_Gamma u = 0,
u(0) = phi_1.
```

The benchmark compares Backward Euler and Crank--Nicolson:

```text
theta = 1.0  Backward Euler
theta = 0.5  Crank--Nicolson
```

## Reproduce

From the project root:

```bash
cmake --build build --target test_heat_graphene_eigenmode_time_convergence
./build/test_heat_graphene_eigenmode_time_convergence
```

The stored output is:

```text
benchmarks/12_heat_graphene_eigenmode_time_convergence/results/time_convergence.csv
```

## Result

The measured errors are:

```text
method,theta,dt,L2_error,rate
BE,1,0.2,0.002067643861895366,
BE,1,0.1,0.001043869101223054,0.9860468996651844
BE,1,0.05,0.0005244883056170111,0.9929583017524573
BE,1,0.025,0.0002628879494996798,0.996462593869439
CN,0.5,0.2,1.103922809190929e-05,
CN,0.5,0.1,2.759508369203772e-06,2.00015613034873
CN,0.5,0.05,6.898532939892645e-07,2.000049768717982
CN,0.5,0.025,1.724567528380127e-07,2.000054966110957
```

This confirms first-order convergence for Backward Euler and second-order
convergence for Crank--Nicolson.
