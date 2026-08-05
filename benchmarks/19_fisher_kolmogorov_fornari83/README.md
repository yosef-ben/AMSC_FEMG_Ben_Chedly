# Benchmark 19: Fisher-Kolmogorov on the 83-Region Connectome

## Goal

Compare the deterministic nodal Fisher-Kolmogorov network model of Fornari et
al. with the P1 metric-graph finite element formulation developed in FEMG.
This benchmark validates the sequential application workflow before any
parallelization or graph reordering is introduced.

## Literature roles

- Weickenmeier et al. provide the one-dimensional sensitivity test used in
  benchmark 18 to validate nonlinear diffusion and reaction.
- Fornari et al. provide the 83-region topology, connectivity-weighted graph
  Laplacian, entorhinal seed, parameters, and biomarker comparison.
- Corti et al. motivate the Fisher-Kolmogorov connectome application and the
  later use of deterministic region-dependent reaction coefficients. The P1
  metric-graph spatial discretization is the contribution of this project.

## Domain verification

The Budapest v3 graph is filtered at five occurrences and aggregated by
FreeSurfer parent region. The resulting topology exactly matches Fornari et
al.: 1015 fine nodes and 37477 fine edges become 83 regions and 1130 edges.
The unweighted degree range is 6-48.

The public Budapest file contains median consensus quantities, whereas
Fornari et al. report averages over 418 subjects. Consequently the topology
matches exactly but the weights are close rather than identical:

```text
quantity                  this benchmark       Fornari et al.
mean connectivity         1.6669               1.57
max connectivity          36.8671              35.32
mean weighted degree      45.3884              42.8
max weighted degree       134.8425             127.6
```

The official mean-mode graph was also downloaded and aggregated with the same
`1%` confidence threshold. It gives mean adjacency `2.4398`, mean weighted
degree `66.4339`, and spectral gap `1.0608`, all farther from the published
values than the median aggregation (`1.6669`, `45.3884`, and `0.7954`). A
unilateral rather than bilateral entorhinal seed does not recover the reported
lobe separation. The audit is stored in
`results/connectome_weight_audit.csv`. This supports retaining the closer
public reconstruction while clearly stating that the original 83-by-83
averaged matrix is unavailable.

## Models

The nodal reference solves

```text
c_dot + L c = alpha c (1-c),
```

where `L=D-A`. The metric-graph FEM solves

```text
M C_dot + H C = R(C).
```

Every coarse connection is represented by one unit P1 element and its
diffusion coefficient is the connectivity weight. Thus `H` has the weighted
graph-Laplacian structure, while `M` and the consistently integrated reaction
term distinguish the FEM model from the nodal reference.

Both models use Backward Euler with Newton, `alpha=0.5`, `dt=0.4` years, and
`T=40` years. The two entorhinal regions start at `c=0.1`; all other regions
start at zero.

## Results

Both solutions remain bounded up to a small transient FEM undershoot of order
`1e-4` caused by the discontinuous nodal seed and the one-element-per-edge
mesh. At `T=40`, both solutions approach the stable state `c=1`.

The biomarker curves have the expected sigmoid shape. Temporal activation is
first, but the four 50-percent crossing times are very close because graph
diffusion rapidly homogenizes the small seed. The ordering of the later lobes
is not yet a quantitative reproduction of Fornari's figure. This is recorded
as a validation limitation, not hidden by parameter tuning.

In `biomarker_comparison`, the upper panels show the four regional curves and
the network average. Since several absolute curves nearly overlap, the lower
panels report each regional curve minus the network average. This makes every
computed trend visible without altering or rescaling the underlying data.

## Discretization study

The spatial study keeps `dt=0.4` and compares 1, 2, 4, and 8 P1 elements per
connection. The maximum difference between the four FEM biomarker curves and
the 8-element reference decreases as follows:

```text
elements per edge    DoFs    maximum difference [%]
1                     83       0.4389
2                   1213       0.1178
4                   3473       0.0229
8                   7993       0.0000 (reference)
```

The small negative transient observed with one element per edge disappears
from two elements per edge onward. This confirms that it is a coarse-mesh
effect and not a failure of the nonlinear solver.

The temporal study uses 8 elements per edge and compares the time steps
`dt=0.8, 0.4, 0.2, 0.1, 0.05`. Successive differences between biomarker curves give the rates

```text
dt       observed rate
0.4      1.0552
0.2      1.0263
0.1      1.0126
```

which converge to the expected first order of Backward Euler. The run with
`dt=0.4` is retained to reproduce the time discretization reported by Fornari
et al.; the refined run with 8 elements per edge and `dt=0.05` is the numerical
reference used to verify the metric-graph FEM implementation.

## Stored output

The `results` directory contains the two biomarker tables, activation times,
the spatial and temporal refinement tables, and report-ready figures.
The generated file `fem_metric_mass.csv` stores the normalized FEM integral,
which is distinct from sampling the solution at the 83 anatomical vertices.
ParaView output remains under `output/fisher_kolmogorov/fornari83/` because it
is generated data.

## Anatomical visualization

The report-oriented output uses each anatomical region exactly once. Open

```text
data/connectome/anatomy/brain_surface.vtk
output/fisher_kolmogorov/fornari83/anatomical/fem_edges.pvd
output/fisher_kolmogorov/fornari83/anatomical/fem_nodes.pvd
```

In ParaView, display the surface in light gray with opacity `0.15-0.25`.
Display `fem_edges.pvd` as thin gray lines. Display `fem_nodes.pvd` using
`Point Gaussian` or a `Glyph` filter with spherical glyphs, color by `c`,
and fix the color range to `[0,1]` for the complete animation. The early
interval `t=0-12` is the informative part of the propagation; at later times
the logistic model approaches its spatially uniform stable state.

The analogous nodal reference files are `nodal_edges.pvd` and
`nodal_nodes.pvd`. Keeping nodes and edges separate avoids duplicated endpoint
values and allows node size and edge thickness to be controlled independently.

For the refined FEM solution, open instead

```text
output/fisher_kolmogorov/fornari83/anatomical_refined/fem_refined_edges.pvd
output/fisher_kolmogorov/fornari83/anatomical_refined/fem_refined_nodes.pvd
```

## Reproduce

Run every command in `commands.txt` from the project root.

## Connectivity diagnostics

The report-ready figure `results/connectome_topology.pdf` reproduces the three
diagnostics used by Fornari et al. and follows the construction of their figure
6: the non-weighted and connectivity-weighted degrees on two stacked glass
brains, and the weighted adjacency matrix, square and as tall as both brains,
with anatomical-group strips along all four of its edges.

Sphere radius and sphere colour encode the same number, that panel's degree,
linearly between that panel's own minimum and maximum. Enlarged nodes therefore
identify network hubs and not anatomical volume. The two panels are normalised
separately, which is the reference's own convention, so a colour in (a) and the
same colour in (b) do not mean the same thing; this is why the degree bar
carries two rows of end numbers rather than one. The connections drawn in the
brains encode nothing: weight appears only in panel (c).

The adjacency uses a base-10 logarithmic colour scale over the full non-zero
range, with nothing clipped at either end. The reason is measurable: the 2260
non-zero entries span 3.603 decades and 77.1% of them fall below 5% of the
maximum, so on a linear ramp three quarters of the connectome would render as
one indistinguishable colour. The consequence is that the top of the range is
compressed, which is why the bar carries decade ticks. Absent connections are
not on the scale at all; the 4629 zero entries, the 83 diagonal ones among
them, are drawn in one flat grey keyed on the bar as "no connection".

Rows and columns keep the solver's own node order, with no permutation,
clustering or thresholding, so that the panel can be laid beside the published
one and compared cell for cell. That is also why the group strips break into 23
runs rather than seven blocks: right hemisphere occupies indices 0 to 40, left
41 to 81, and the brainstem 82.

For the public reconstruction used by the solver:

```text
unweighted degree range = 6 to 48
weighted degree range   = 2.2885 to 134.8425
adjacency range         = 0.0092 to 36.8671
```

The ranges and block structure are qualitatively consistent with the
small-world organization reported in the paper. Exact weighted values are not
expected to coincide because the original matrix is an average over 418
subjects, whereas this benchmark uses the downloadable median consensus
graph. Per-region values are stored in
`results/connectivity_diagnostics.csv`.

## Literature verification

Checked against Fornari, Schafer, Jucker, Goriely and Kuhl, J. R. Soc.
Interface 16 (2019) 20190356, sections 3 and 4.1.

Confirmed against the paper:

| Quantity | Fornari et al. | This benchmark |
|---|---|---|
| fine graph | 1015 nodes, 37477 edges | 1015, 37477 |
| coarse graph | 83 nodes, 1130 edges | 83, 1130 |
| unweighted degree range | 6 to 48 | 6 to 48 |
| mean fibre number | 40.2 | 40.1619 |
| mean fibre length | 38.40 mm | 38.4009 mm |
| mean adjacency | 1.57 | 1.6669 |
| mean weighted degree | 42.8 | 45.3884 |
| adjacency range | 0.01 to 35.32 | 0.0092 to 36.8671 |
| weighted degree range | 2.1 to 127.6 | 2.2885 to 134.8425 |
| model | `dc/dt = -L c + alpha c (1-c)`, equation (3.2) | same |
| conversion rate | `alpha = 0.5` | 0.5 |
| seed | entorhinal cortex at `c_0 = 0.1` | both entorhinal vertices |
| time integration | implicit, 100 steps of `dt = 0.4` years | same |
| biomarker | `C(t) = sum_I c_I(t)`, equation (4.1) | same, normalised to a percentage |

The paper defines the edge weight as the mean fibre number divided by the mean
fibre length, `A_IJ = n_IJ / l_IJ`. The preprocessing here instead sums the
Budapest `electrical_connectivity_median` of the fine connections. Building the
Laplacian both ways gives

```text
weighting                    mean adjacency   mean weighted degree   lambda_2
electrical (this benchmark)         1.6669                45.3884     0.7954
n_IJ / l_IJ (paper formula)         1.3978                38.0599     0.7431
paper reported values               1.57                  42.8            --
```

so the current choice is the closer of the two to the published statistics and
is retained.

### The one quantity that is not reproduced

Figure 7 of the paper shows the four lobe biomarkers separated by roughly five
years, activating in the order temporal, frontal, parietal, occipital. Here
they cross the 50-percent level within `2.9e-7` years of each other.

Benchmark 23 resolves this. With `lambda_2 = 0.795` and `alpha = 0.5` the
Damkohler number `alpha / lambda_2` is `0.63`: the connectome homogenises
faster than the reaction grows, so coincident lobe curves are the correct
solution of equation (3.2) at unit weight scale. Reproducing the published
separation requires scaling the Laplacian by a factor near `0.005`. The
published figure therefore cannot come from the connectivity-weighted
Laplacian used at unit scale, whatever the reconstruction of the weights.

The activation order is also worth recording. At any scaling that separates the
lobes, the frontal lobe is last here, because the frontal pole is the most
weakly connected node of the graph. The paper reports the same property in its
figure 9 discussion, where the frontal pole has the longest infection time of
all 83 seeding regions, so our ordering is consistent with the paper's own
network diagnostics even though it differs from its figure 7.

### Time discretization

The retained run uses `dt = 0.4` to match the paper. The temporal refinement
table in this benchmark shows the global 50-percent crossing moving from
`11.50` years at `dt = 0.8` to `13.70` years at `dt = 0.05`, so the headline
value carries a time-discretization error of order one year. The paper uses the
same step with a first-order implicit scheme and therefore carries a comparable
error; the agreement of the absolute crossing times should be read with that in
mind.

All CSV outputs in `results/` were regenerated from source in this working copy
and reproduce byte for byte.
