# MPI condensation

The condensed algorithm of benchmark 29 distributed with MPI. The
connectome is not vertex-partitioned: every rank owns one contiguous block
of edges, balanced by interior unknowns, and with it their interior DoFs;
the 83 original vertices are the replicated interface. Per step every rank
condenses its own edges into a partial vertex system, two MPI_Allreduce
calls sum the 83x83 Schur complement and the condensed right-hand side (a
constant payload of about 55 kB, independent of the refinement), every
rank solves the interface system redundantly and back-substitutes its own
interiors. MPI is enabled for the one new executable alone
(`test_condensation_mpi`).

## Validation

`results/mpi_validation.csv`: the distributed state, gathered at every
step, against the sequential condensed engine advanced on rank 0:

```text
cells   8, ranks 2:  max difference 7.2e-16
cells   8, ranks 4:  max difference 7.5e-16
cells   8, ranks 8:  max difference 8.9e-16
cells 128, ranks 8:  max difference 5.5e-15
```

## Strong scaling

128 cells per edge (143593 unknowns), 200 steps, three processes per rank
count with per-field medians. Up to four ranks each process is bound to
its own physical core; eight ranks use the hardware threads of the 4-core,
8-thread laptop CPU, matching the OpenMP study. Phase medians are the
maximum across ranks; no barrier separates the collective from the wait it
absorbs, so the Allreduce column carries both the wire time and the
per-step imbalance. The sequential reference T1 is pooled from repeated
runs of the sequential engine in the same session, 8.07 ms per step.

```text
ranks               step [ms]   local  allreduce  iface   back |     S       E
sequential baseline      8.07                                  |  1.00      -
      2                  4.66     4.30       0.09   0.07   0.21 |  1.73    0.87
      4                  4.29     3.24       0.71   0.11   0.23 |  1.88    0.47
      8                  2.56     1.98       0.20   0.14   0.25 |  3.15    0.39
```

The MPI implementation executed with one rank takes 8.93 ms per step, 11
percent above the sequential engine; its row is kept in the CSV and
speedups are computed against the sequential engine.

## Reading, next to the OpenMP record

```text
              2 workers   4 workers   8 workers
OpenMP step    5.65 ms     3.37 ms     2.06 ms
MPI step       4.66 ms     4.29 ms     2.56 ms
```

1. On this single machine shared memory wins: OpenMP merges the
   thread-local vertex systems in 0.01-0.05 ms, while the MPI collective
   costs 0.09-0.71 ms per step and, with all physical cores busy, absorbs
   the per-step imbalance of the fully loaded laptop (the 4-rank column).
   Communicating at every one of the 200 steps, the distributed version
   pays that price 200 times.
2. The local edge work itself scales the same way in the two models
   (3.24 ms on 4 ranks against 3.06 ms on 4 threads): the algorithm is
   insensitive to how the workers are implemented, and the difference
   lives entirely in the synchronization.
3. The MPI structure is nevertheless the one that leaves the machine:
   the payload of the only collective is constant in the refinement, the
   interface solve is redundant instead of communicated, and the
   partition owns its interiors outright, so the same code addresses
   multiple nodes where shared memory cannot go. On one node, OpenMP is
   the appropriate choice.

## Files

- `results/mpi_validation.csv` - distributed against sequential, every step.
- `results/mpi_scaling.csv` - phases, pooled T1, speedup and efficiency.
