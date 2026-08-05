# Linear Solver Complexity Note

This additional benchmark compares Eigen linear solvers on the elliptic matrix

```text
A = H + M,
```

corresponding to the model problem `-Delta_Gamma u + u = f` on the
graphene-like graph.

The tested methods are:

```text
CG
GMRES
BiCGSTAB
SparseLU
```

Eigen 3.3.9 does not provide a QMR implementation. For this reason the report
should not include a curve labelled `QMR` unless a dedicated implementation is
added. `BiCGSTAB` is included as the available nonsymmetric Krylov comparison.

The reference curves in `linear_solver_complexity.svg` are:

```text
N_dof
N_dof^2
N_dof^3
```

The benchmark data are stored in:

```text
results/graphene_linear_solver_complexity.csv
```

The figure is stored in:

```text
results/linear_solver_complexity.svg
```
