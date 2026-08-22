# Benchmark suite

Each directory is one self-contained experiment: a `README.md` stating the
problem, the parameters and the reference it follows, a `commands.txt` that
regenerates it from scratch, and the report-ready results.

`FIGURES.md` records, for every figure used in the report, the script that
draws it, the file it reads and the executable that wrote that file.
`scripts/verify-figures.py` re-derives from the stored results the quantities
the figures and these READMEs state, and fails if any of them disagree.

## Order

The suite is meant to be read in order: each group depends on the one before
it.

### Linear validation, 02 to 16

The four-pointed star of Abbate and Di Primio with constant, linear and sine
solutions, a genuinely time-dependent radial decay, an eigenmode relaxation,
then spectral problems and heat equations on the graphene-like and binary-tree
graphs, with spatial and temporal convergence, a spectral comparison, a
complexity study and an energy-decay check.

### Nonlinear solver, 17 and 18

`18` is the validation that matters: the one-dimensional sensitivity test of
Weickenmeier et al., with the front-speed check against the Fisher-KPP
asymptotic speed and the time-step study. `17` predates it and is retained for
reproducibility only; it is superseded by `21`.

### Connectome application, 19 to 21

`19` reconstructs and audits the 83-region graph of Fornari et al. and compares
the nodal network model with the metric-graph formulation, with spatial and
temporal refinement. `20` is the conversion-rate study. `21` is the
deterministic model of Corti et al. with region-dependent reaction rates.

### Analysis, 23

`23` establishes the dimensionless group that governs every connectome result,
`Da = alpha / (rho lambda_2)`. It explains why the lobe biomarkers of `19`
coincide while those of `21` separate, and it bounds the consistent-mass P1
discretization, which loses boundedness above `Da` of about 13. Read it
immediately after `19`, since it answers the question `19` leaves open.

### Domain and dynamics, 24 to 27

`24` documents the graph itself, `25` the seeding study (a record: in the
finite element model the ranking follows the seed mass, see its README),
`26` the anatomical progression and the activation-time map, and `27` the
staged spreading of the two clinical seedings, tau and amyloid-beta, against
the progression the literature expects.

### Performance, 22

The optimised single-process baseline, measured on the problem of `21`. It is
the starting point of the reordering and parallel work and should be read last.

## Where to start

To retrace the whole chain from scratch:

```bash
python3 scripts/prepare-fornari-connectome.py     # rebuild the graph
python3 scripts/audit-fornari-connectome.py       # audit it against the paper
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
./build-release/test_fisher_kolmogorov_logistic   # orders: semi-implicit
                                                  # 1.99918, 1.99966;
                                                  # fully implicit
                                                  # 1.00144, 1.00036
/usr/bin/python3 scripts/verify-figures.py        # every figure against its data
```

Then run `commands.txt` in `18`, `19`, `23`, `20`, `21`, `25`, `26`, `22`, in
that order.
