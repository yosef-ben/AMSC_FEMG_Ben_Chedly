# Benchmark 22: Sequential Fisher-Kolmogorov Performance

## Goal

Establish the optimized single-process baseline immediately before graph
reordering and parallelization. The benchmark solves exactly the deterministic
Corti-83 problem from benchmark 21 while refining every metric connection.

Visualization, CSV callbacks, and VTK output are disabled during timing. The
scientific model, coefficients, initial condition, time step, and number of
time steps are unchanged.

## Experimental protocol

The executable is compiled in a separate Release build:

```text
compiler       GNU C++ 13.1.0
build type     Release (-O3 -DNDEBUG)
CPU            Intel Core i5-10210U, 4 cores / 8 hardware threads
threads used   1
T              20 years
dt             0.2 years
time steps     100
warm-up        1 run per mesh
measurements   median of 5 runs per mesh
```

`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS` are fixed
to one by the runner. The raw measurements are retained, since laptop
frequency scaling and operating-system activity produce visible run-to-run
variation.

## Measured phases

- **Initialization:** graph input and global DoF numbering.
- **Assembly:** interpolation of coefficients and construction of mass and
  diffusion matrices.
- **Solve:** 100 semi-implicit Crank-Nicolson steps.
- **Total:** input, setup, assembly, and solve, excluding process startup.

At every nonlinear step, the extrapolated saturation coefficient changes.
Consequently the sparse time-step matrix is rebuilt and factorized at every
step. The current `Eigen::SparseLU` type uses Eigen's default
`COLAMDOrdering<int>`; this is the baseline ordering against which future RCM
experiments must be compared.

## Results

| Cells/edge | DoFs | Matrix nnz | Assembly [s] | 100 steps [s] | Time/step [s] |
|---:|---:|---:|---:|---:|---:|
| 1 | 83 | 2343 | 0.000563 | 0.091428 | 0.000914 |
| 2 | 1213 | 5733 | 0.001679 | 1.739028 | 0.017390 |
| 4 | 3473 | 12513 | 0.002690 | 1.738110 | 0.017381 |
| 8 | 7993 | 26073 | 0.004135 | 2.404785 | 0.024048 |

Matrix nonzeros and assembly cost grow regularly with refinement. By contrast,
the direct-solver time is not a function of the DoF count alone. The
one-interior-node mesh has a hub-dominated bipartite sparsity pattern and is
not cheaper than the four-cell mesh despite having fewer unknowns. Sparse LU
factorization depends on elimination ordering and fill-in as well as matrix
dimension.

For the refined cases, more than 99% of the measured total time is spent in
the 100 time steps; assembly is negligible. Therefore the first optimization
target is the repeated sparse factorization, not finite element assembly or
file input.

These four points are a baseline study, not sufficient evidence for a single
asymptotic complexity law. The non-monotone result is retained rather than
smoothed because it motivates measuring bandwidth, profile, factor fill, and
factorization time in the reordering study.

## Stored output

`results/sequential_performance.csv` contains medians and min/max intervals.
`results/raw_timings.csv` contains all 20 measured runs.
`results/sequential_performance.png` and `.pdf` are report-ready figures.

## Reproduce

From the project root:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target test_fisher_kolmogorov_corti83
python3 scripts/run-fisher-sequential-performance.py \
  --executable build-release/test_fisher_kolmogorov_corti83 \
  --output-dir output/fisher_kolmogorov/corti83/performance
python3 scripts/plot-fisher-sequential-performance.py \
  output/fisher_kolmogorov/corti83/performance/sequential_performance.csv \
  --output output/fisher_kolmogorov/corti83/performance/sequential_performance.png
```

## Next boundary

This benchmark deliberately stops before optimization. The next stage should
measure, for identity, COLAMD, and RCM orderings:

1. matrix bandwidth and profile;
2. symbolic and numerical factorization time;
3. nonzeros in the LU factors;
4. end-to-end time for the 100-step simulation.

Only after the best sequential ordering is established should strong and weak
parallel-scaling results be interpreted.

## Literature context

This benchmark has no counterpart in the reference papers, which report timings
only in passing. One comparison is still available: Fornari et al. state that
their 83-node nodal Fisher-Kolmogorov simulation "runs 0.55 and 0.64 s without
and with output on a standard laptop computer" for 100 implicit steps of
`dt = 0.4` years.

The corresponding metric-graph FEM run here, at one element per edge and the
same 83 unknowns, takes `0.091` s for its 100 semi-implicit steps. The two
machines are different and the schemes are not identical, so this is an
order-of-magnitude check rather than a benchmark, but it confirms that the P1
metric-graph formulation carries no cost penalty over the nodal network model
at the coarsest mesh, while allowing arbitrary refinement inside the edges.

The degrees of freedom and matrix nonzeros in the results table were
regenerated from source in this working copy and reproduce exactly. The
timings are machine-dependent and were not re-measured.

## Figure note

The figure no longer carries `O(N)` and `O(N^2)` guide lines. Four points on a
sequence of meshes whose sparsity pattern changes with refinement do not
identify a complexity law, and the guides contradicted the text above. They are
replaced by the measured minimum-to-maximum band over the five repetitions,
which is the quantity the protocol actually establishes.
