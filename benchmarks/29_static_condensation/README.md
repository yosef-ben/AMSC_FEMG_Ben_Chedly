# Static condensation of the semi-implicit step

The ordering study (benchmark 28) showed that the interiors-first DoF
numbering of the library is implicit static condensation: eliminating the
edge chains first leaves fill only in the vertex block. This record makes
that elimination explicit. The semi-implicit time step

```text
[M + (dt/2) H - (dt/2) W(c_hat)] c_new = [M - (dt/2) H + (dt/2) W(c_hat)] c_old
```

is advanced without ever assembling a global matrix: every cell contributes
a 2x2 block (consistent mass, the 3-point Gauss reaction weight of the
library), the blocks are accumulated per edge into its interior tridiagonal
system, a Thomas elimination condenses each edge to a 2x2 contribution to
the Schur complement on the 83 original vertices, the dense 83-vertex
system is solved once, and the interiors are recovered by one saxpy per
edge from the vectors stored during the elimination. The per-edge loop
touches only edge-local buffers and the vertex accumulators: it is the
parallel unit of the following stages.

`test_condensation` (test_problem/performance, with the stepper in
`condensed_stepper.hpp`) runs the deterministic Corti-83 reference problem
of the sequential-performance record; the production solver is not
modified.

## Validation

`results/validation.csv`: the condensed trajectory against the `solve()`
of the production class, at every one of the 100 steps:

```text
cells  1: max difference 6.2e-15
cells  2: max difference 1.0e-14
cells  4: max difference 9.7e-15
cells  8: max difference 4.2e-15
```

The requirement was 1e-10; the agreement is at round-off, and since the
class run reproduces the stored benchmark records byte for byte (their own
guards), the condensed solver reproduces the stored biomarker and crossing
results transitively.

## Results

Median per-step times (`results/condensation_study.csv`; full = the
SimplicialLDLT full-system loop on the natural ordering, the reference
selected by benchmark 28, measured in the same process):

```text
cells      n   condensed   local  iface   back |   full   ratio
    8   7993     1.03 ms    0.90   0.11   0.02 |   3.61    3.5x
   16  17033     1.99 ms    1.84   0.11   0.04 |   8.34    4.2x
   32  35113     4.11 ms    3.89   0.12   0.10 |  13.52    3.3x
   64  71273     6.73 ms    6.45   0.10   0.18 |  35.02    5.2x
  128 143593    10.99 ms   10.56   0.09   0.34 |  64.45    5.9x
```

The final states of the two engines agree within 2.8e-14 at every size.

## Reading

1. The explicit condensation beats the best full-system configuration by
   3.5x to 5.9x, growing with the refinement, because it removes the
   sparse-matrix machinery (triplets, compressed storage, symbolic
   analysis) that the ordering study identified as the residual cost.
2. The interface solve is 0.1 ms and does not grow with the refinement:
   the serial fraction of the step falls from 11 percent at 8 cells to
   0.8 percent at 128, so the per-edge loop, 96 percent of the step at
   128 cells, is the parallel content of the following stages, with an
   Amdahl ceiling of about 7.5 on 8 cores at that size.
3. Cumulative algorithmic gain at 64 cells per edge, to be reported
   separately from any parallel speedup: production SparseLU+COLAMD
   140.3 ms per step (benchmark 28), full-system LDLT on the natural
   ordering 21-35 ms, condensed 6.7 ms.
4. Scaling size: 128 cells per edge (143593 unknowns) is the smallest
   refinement with a sequential step above 10 ms, the practical floor
   for measuring 8-way strong scaling; the topology is kept fixed and
   only the finite-element resolution along each connection is
   increased, providing a controlled increase in the number of unknowns
   without modifying the underlying connectome.

## Files

- `results/validation.csv` - condensed against the production class.
- `results/condensation_study.csv` - phases of both engines and their
  final-state difference.
