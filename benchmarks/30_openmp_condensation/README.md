# OpenMP condensation

The per-edge condensation of benchmark 29 parallelized with OpenMP. The
mathematics is unchanged: every thread condenses a static block of edges
into its own copy of the 83-vertex system, one critical section per thread
merges the copies, the dense interface solve stays serial and the
back-substitution loop is parallel again. The work vectors are indexed by
the edge offsets, so distinct edges never share a slice and the loops are
race-free by construction; only the summation order of the vertex system
depends on the thread count.

OpenMP is enabled for the one new executable alone
(`test_condensation_omp`); the library and every other target build exactly
as before.

## Validation

`results/omp_validation.csv`: the OpenMP engine against the sequential
condensed engine, advanced in lockstep and compared at every one of the
200 steps:

```text
cells   8, threads 2:  max difference 2.3e-15
cells   8, threads 4:  max difference 2.8e-15
cells   8, threads 8:  max difference 2.7e-15
cells 128, threads 8:  max difference 2.5e-14
```

## Strong scaling

128 cells per edge (143593 unknowns), 200 steps per run, three processes
per thread count with per-field medians, threads pinned
(`OMP_PLACES=cores`, `OMP_PROC_BIND=close`). The machine is a 4-core,
8-thread laptop CPU (Intel i5-10210U) whose turbo makes single-process
timings drift by up to 30 percent, so the sequential reference `T1` is the
pooled median of the sequential runs of every process, 9.58 ms per step
(`results/omp_scaling.csv` keeps every raw field).

```text
threads               step [ms]   local  reduce  iface   back |     S       E
sequential baseline        9.58                                |  1.00      -
      2                    5.65     5.31    0.01   0.08   0.26 |  1.69    0.85
      4                    3.37     3.06    0.02   0.08   0.20 |  2.84    0.71
      8                    2.06     1.72    0.05   0.09   0.21 |  4.64    0.58
```

The OpenMP implementation executed with a single thread takes 10.70 ms per
step, 12 percent more than the sequential engine (thread-local
accumulators and barriers); its row is kept in the CSV and speedups are
computed against the sequential engine, never against the one-thread
OpenMP run and never against the retired production solver.

## Reading

1. The clean strong-scaling result is the one on the four physical cores:
   S = 2.84 at 71 percent efficiency for the whole step, with the local
   edge phase alone scaling from 10.25 to 3.06 ms (84 percent parallel
   efficiency on the parallelizable part; the two figures refer to
   different scopes and are both correct).
2. Eight threads are four cores with simultaneous multithreading, not
   eight physical cores: the step still falls to 2.06 ms (S = 4.64), and
   the reduced efficiency at eight threads is consistent with the use of
   SMT on a four-core processor and with the increasing relative cost of
   the serial and memory-bound phases.
3. Those parts are now visible and small: the interface solve (0.08 ms,
   constant by design), the reduction of the thread-local vertex systems
   (0.05 ms at 8 threads) and the back-substitution, a pure memory stream
   that saturates the shared bandwidth near 0.20 ms.
```text
algorithmic gain (benchmark 28-29):  production 64.5-140.3 ms -> condensed 9.6 ms
parallel gain (this record):         condensed 9.6 ms -> 2.06 ms on 8 threads
```

## Files

- `results/omp_validation.csv` - lockstep OpenMP against sequential.
- `results/omp_scaling.csv` - phases, pooled T1, speedup and efficiency.
