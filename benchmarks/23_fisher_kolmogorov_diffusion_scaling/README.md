# Benchmark 23: Diffusion Scaling of the 83-Region Connectome

## Goal

Resolve the open validation limitation of benchmark 19. That benchmark
reproduces the topology, the parameters and the time discretization of Fornari
et al. exactly, yet its four lobe biomarker curves activate within `1e-7`
years of each other, whereas figure 7 of the paper shows a separation of
roughly five years between the temporal and the occipital lobe.

This benchmark shows that the discrepancy is neither a defect of the
implementation nor a property of the reconstructed connectome. It is entirely
controlled by one modelling quantity, the scale of the connectivity weights,
and it is governed by a single dimensionless number.

## The dimensionless number

The nodal model of Fornari et al. is

```text
dc/dt = -rho L c + alpha c (1-c),
```

where `L` is the connectivity-weighted graph Laplacian and `rho` is a uniform
scaling of the weights. Writing `rho = 1` recovers equation (3.2) of the paper
literally.

Linearising around a small uniform concentration, the mean over the graph grows
at the reaction rate `alpha`, while any lobe-to-lobe imbalance relaxes at the
rate `rho * lambda_2`, with `lambda_2` the Fiedler value of `L`. Regional
separation therefore survives only when the reaction is faster than graph
homogenisation, that is when the Damkohler number

```text
Da = alpha / (rho * lambda_2)
```

is large. For the reconstructed Budapest-83 graph:

```text
lambda_2      = 0.772254
lambda_max    = 145.4535
mean weighted degree = 42.7547
max adjacency w      = 35.3221
```

## The same number for the metric-graph FEM

`Da` above is the Damkohler number of the nodal model, whose mass matrix is
the identity. The finite element model shares `L` exactly at one element per
connection, but every vertex also carries half of each incident connection:
with unit lengths the lumped mass is `M_ii = deg_i / 2`, from 3 to 24 and
13.61 on average, and the consistent mass has the blocks h/3, h/6. The rate
of the slowest transport mode is then the smallest non-zero generalized
eigenvalue of `(L, M)` instead of the Fiedler value of `L`
(`scripts/tabulate-connectome-mass-spectrum.py`, `results/mass_spectrum.csv`):

```text
mass matrix   lambda_2    ratio to the identity   Da at rho = 1 / 0.05 / 0.005
identity      0.772254    1.00                    0.65 /  12.9 /  129
lumped        0.060627   12.74                    8.25 / 164.9 / 1649
consistent    0.063533   12.16                    7.87 / 157.4 / 1574
```

At the same `rho` the finite element model therefore sits at a Damkohler
number 12 to 13 times larger than the nodal one, which is why its lobes
separate at a smaller nominal `Da` in the report figure: the 6.2 years of the
consistent mass at nominal `Da = 12.9` correspond to `Da` of about 160 for
that model, where the nodal sweep below gives between 4.38 (Da 129) and 5.90
(Da 324) years. The two models agree on what controls the separation; the
number at which a given separation is reached depends on the mass matrix,
as the 1.5-year shift of benchmark 19 already does.

## Measured behaviour

The nodal reference of `test_fisher_kolmogorov_fornari83` was run for
`alpha = 0.5`, `dt = 0.4` years and `T = 120` years, sweeping `rho`. The lobe
spread is the difference between the latest and the earliest 50-percent
crossing time among the temporal, frontal, parietal and occipital lobes.

| rho | Da | lobe spread [yr] | metric-graph FEM |
|---:|---:|---:|:---|
| 1.0 | 0.65 | 3.87e-07 | bounded |
| 0.5 | 1.29 | 2.84e-05 | bounded |
| 0.2 | 3.24 | 0.0080 | bounded |
| 0.1 | 6.47 | 0.1317 | bounded |
| 0.05 | 12.95 | 0.6529 | bounded |
| 0.04 | 16.19 | 0.9210 | leaves [0,1] |
| 0.03 | 21.58 | 1.3204 | leaves [0,1] |
| 0.025 | 25.90 | 1.5986 | leaves [0,1] |
| 0.02 | 32.37 | 1.9586 | leaves [0,1] |
| 0.01 | 64.75 | 3.1600 | leaves [0,1] |
| 0.005 | 129.49 | 4.3801 | leaves [0,1] |
| 0.002 | 323.73 | 5.8966 | leaves [0,1] |
| 0.001 | 647.46 | 6.9863 | leaves [0,1] |

The spread varies over seven decades while the topology, the seed and the
reaction rate are unchanged. This is the quantitative statement that benchmark
19 was missing.

## Consequences for the earlier benchmarks

- **Benchmark 19** uses the literal `rho = 1`, giving `Da = 0.65`. The
  connectome homogenises faster than the reaction grows, so the four lobe
  curves must coincide. Its near-zero lobe separation is the correct solution
  of the equation as we reconstruct it, not a validation failure. A
  separation of the magnitude of figure 7 of the paper is obtained at `Da` of
  order `10^2`, that is `rho` near `0.005`, about two hundred times weaker
  than the literal scale. The paper does not state the scale of its Laplacian
  relative to the time unit (the adjacency is fibres per millimetre, alpha
  and dt are per year), and neither its 83x83 matrix nor its code is
  distributed; so this benchmark establishes what the separation depends on,
  not which choice the paper made. See the report's "Limits of an exact
  reproduction" for the full list of unspecified details.
- **Benchmark 21** normalises the diffusivity as `D_e = w_e/max(w)`, that is
  `rho = 1/35.3221 = 0.028311`, and uses the seven Corti reaction means whose
  average over the 83 vertices is `0.1252`. Its Damkohler number is therefore
  `5.73`, in the regime where regional curves separate. The clearly ordered
  regional averages of benchmark 21 and the coincident lobe curves of
  benchmark 19 are consistent with each other: the two benchmarks sit on
  opposite sides of the same crossover.

## Validity boundary of the metric-graph FEM

The sweep also runs the P1 metric-graph FEM at one element per edge, with both
available time schemes, and classifies a run as bounded when its final state
lies in `[0,1]`. The violation of the physical range grows with `Da` well
before the failure: the transient undershoot is `-1.3e-4` at `Da = 0.65`
(the stored biomarkers of benchmark 19) and already `[-0.17, 1.04]` at
`Da = 12.9` (`results/fem_transient_rho_0p05.csv`), from which the solution
still recovers. The solution is bounded up to `rho = 0.05` (`Da = 12.9`) and
unbounded from `rho = 0.04` (`Da = 16.2`) onwards, the three scalings 0.04,
0.03 and 0.025 having been added to locate the boundary; beyond it the
solution no longer recovers and diverges:

```text
rho = 0.005, one element per edge, Corti scheme:
  nodal reference  [0.999848, 0.999997]
  metric-graph FEM [-506.980,  523.169]
```

The boundary therefore separates recoverable from unrecoverable violations
of the maximum principle, not violation from compliance. The nodal reference
stays in `[0,1]` at every scaling tested. The failure occurs
identically with Backward Euler and Newton and with the semi-implicit
Crank-Nicolson scheme, so it is not a nonlinear-solver failure. It also
persists at 2, 4 and 8 elements per edge, so it is not a matter of resolving
the front.

The cause is the consistent P1 mass matrix. Its off-diagonal entries are
positive, so the semi-discrete system does not satisfy a discrete maximum
principle: the reaction at one node, multiplied by the inverse of the mass
matrix, can drive a neighbouring node below zero. While diffusion dominates,
the diffusion matrix, an M-matrix, compensates; once the reaction dominates
nothing does, and Newton then converges to the wrong branch of the logistic
term. The standard remedy is a diagonal mass, and it is available in
`fisher_kolmogorov_problem` through `set_mass_lumping(true)`: the row-sum
lumped mass (`h/2` per node and cell) with the reaction and the reaction
weight matrix evaluated by the vertex rule, so that the reaction at a node
acts on that node alone. The nodal model of Fornari et al. and the Hadamard
product in equation (4) of Corti et al. are of the same kind, with the
identity as mass; the lumped mass differs from it by the per-vertex factor
`deg/2` at unit edge lengths, a rescaling of time vertex by vertex.

## The lumped variant over the sweep

`results/fem_lumped_sweep.csv`, written by
`scripts/tabulate-fem-lumped-sweep.py` from the stored biomarker files,
records the mass-lumped FEM at every scaling of the sweep with both schemes,
plus the four-element case at `rho = 0.005` where the consistent mass diverged
(`test_fisher_kolmogorov_fornari83 ... be_lumped | cn_lumped`):

```text
rho     Da      lumped BE: range over time   spread   lumped CN: range   spread
1.0     0.65    [0, 1]                       0.027    [0, 1]             0.028
0.05    12.95   [0, 1]                       5.09     [0, 1]             5.86
0.04    16.19   [0, 1]                       5.54     [0, 1]             6.36
0.005   129.5   [0, 1]                       9.37     [0, 1]             10.58
0.001   647.5   [0, 1]                       12.11    [0, 1]             13.62
0.005, four elements per edge: [0, 1] with both schemes (consistent mass: [-507, 523])
```

Every one of the 28 runs stays within `[0, 1]` to machine precision over the
whole transient. At the literal scale the lumped and the consistent mass
agree, network crossing 12.72 against 12.68 years, and both keep the four
lobes together (spreads 0.027 and 0.008 years); at large `Da` the lumped
formulation separates the lobes more than the nodal model does, because its
mass `deg/2` slows the diffusion at the well-connected vertices relative to
the identity mass, a difference between two diagonal-mass discretizations
and not a defect of either. The comparisons of the report keep the
consistent mass, the formulation of the first part of the work, which is
valid at every `Da <= 5.73` used there; the lumped option is the extension
into the reaction-dominated regime.

Benchmarks 19 and 21 are unaffected: both sit at `Da <= 6`, well inside the
bounded region.

## Reproduce

Run the commands in `commands.txt` from the project root.

## Stored output

```text
results/diffusion_scaling.csv           sweep table, including FEM status
results/diffusion_scaling_summary.json  spectral quantities and reference scalings
results/diffusion_scaling.png/.pdf        report figure: nodal, FEM consistent mass, FEM lumped mass at three scalings
results/diffusion_scaling_summary_rows.csv  extremes, stopping time and lobe spread of its nine runs
results/diffusion_scaling_spread.png/.pdf   record figure: nodal biomarkers plus spread against Da over the sweep
results/fem_lumped_sweep.csv              the lumped FEM over the whole sweep
results/mass_spectrum.csv                 Fiedler value of L against the identity, lumped and consistent mass
```

The report figure shows, under identical conditions (one element per
connection, fully implicit scheme, dt = 0.4, T = 40, entorhinal seed), the
nodal model, the consistent-mass FEM and the lumped-mass FEM at rho = 1, 0.05
and 0.005, with the envelope of the 83 vertex concentrations shaded. The
lobes separate as Da grows in all three: spreads 3.9e-7, 0.65 and 4.4 years
(nodal), 0.008 and 6.2 years (consistent, then breakdown), 0.027, 5.1 and 9.4
years (lumped). With the consistent mass the envelope leaves [0,1] at
Da = 12.9 (-0.08, 1.02) and the Newton iteration fails at t = 15.6 years at
Da = 129 (-0.41, 1.17), where the run stops; with the lumped mass every panel
stays within [0,1].

The record figure `diffusion_scaling_spread` shows the nodal biomarkers at
three scalings and the measured spread against the Damkohler number, with the
setting of benchmark 19 and the separation reported by Fornari et al. marked
and the scalings at which the consistent-mass FEM leaves [0,1] shaded; the
report quotes the spread values in the text instead. The
lobe colours follow figure 7 of the paper; every pair of them separates by at
least `10.8` in OKLab under protan, deutan and tritan simulation, and the line
styles repeat the identity so it never rests on colour alone.

## What this benchmark does not claim

The activation order at large `Da` is temporal, then occipital, then parietal,
then frontal, whereas Fornari et al. report temporal, frontal, parietal,
occipital. The frontal lobe is last here because the frontal pole is the most
weakly connected node of the graph, which the paper itself reports in its
figure 9 discussion of infection times. Matching the published order would
require a different lobe assignment or a different weighting, and is not
attempted.
