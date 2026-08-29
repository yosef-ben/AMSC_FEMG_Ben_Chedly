# Hybrid MPI+OpenMP comparison

One question, kept deliberately small: at equal total worker count, does
combining distributed and shared memory help on this machine? The hybrid
binary (the MPI driver with the per-rank edge loops under OpenMP) runs the
same eight hardware threads of the 4-core machine in four arrangements,
in one session, against the pooled sequential reference of 8.40 ms per
step, on the 128-cell mesh (143593 unknowns, 200 steps, three repeats,
ranks packed with --map-by slot:PE=threads and threads bound inside).

## Validation

`results/hybrid_validation.csv`: lockstep against the sequential engine,
every step: 7.8e-16 (8 cells, 2x4), 8.9e-16 (8 cells, 4x2), 5.1e-15
(128 cells, 2x4).

## Comparison

```text
ranks x threads   step [ms]   local  allreduce  iface   back |     S       E
    1 x 8            2.06      1.78       0.00   0.09   0.19 |  4.08    0.51
    2 x 4            2.30      1.91       0.07   0.09   0.23 |  3.65    0.46
    4 x 2            2.35      1.88       0.13   0.11   0.23 |  3.57    0.45
    8 x 1            2.61      2.01       0.21   0.13   0.25 |  3.22    0.40
```

The 1x8 and 8x1 rows reproduce the pure OpenMP (2.06 ms, benchmark 30)
and pure MPI (2.56 ms, benchmark 31) measurements of their own records.

## Reading

The local edge work is the same in every arrangement (1.78 to 2.01 ms):
the ordering is decided entirely by the synchronization, which grows
monotonically with the rank count (Allreduce 0.00 to 0.21 ms, redundant
interface solve 0.09 to 0.13 ms). On one machine, every rank added at
fixed worker count converts a free shared-memory merge into paid
communication, so pure OpenMP is the best arrangement and the hybrids
interpolate between the two purities. Combining the two models offers no
advantage on a single node; its role is the multi-node setting the
machine cannot exercise, threads inside a node and ranks across nodes,
for which this binary is the working implementation.

## Files

- `results/hybrid_validation.csv` - lockstep against sequential.
- `results/hybrid_comparison.csv` - the four arrangements, phases, S and E.
