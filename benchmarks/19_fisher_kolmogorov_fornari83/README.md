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

Two provenance facts are not stated in the paper and are recorded here so the
reconstruction is honest about what is choice and what is data. First, the
five-occurrence threshold is not declared by Fornari et al.; it is the only
integer that reproduces their fine graph of 1015 nodes and 37477 edges, since
four keeps 40895 edges and six keeps 34718. Second, although the paper
describes the fibre quantities as cohort means over 418 subjects, the
published graph statistics are reproduced by the median fibre-number and
fibre-length fields of the public dataset, not by its mean fields. Weighting
every fine connection by `fiber_count_median / fiber_length_median` and
summing parallel connections, as conductances in parallel, gives

```text
quantity                  this benchmark       Fornari et al.
mean connectivity         1.5702               1.57
max connectivity          35.3221              35.32
mean weighted degree      42.7547              42.8
max weighted degree       127.6435             127.6
min weighted degree       2.0505               2.1
```

every value at the precision at which the paper prints it. The regions at the
extremes match as well: weighted degree from the right frontal pole to the
right precentral gyrus, largest adjacency and largest fibre count both on the
right superior parietal to precuneus pair, smallest adjacency on the left
lateral orbitofrontal to isthmus cingulate pair. The paper names the same
regions, although its sentence attaches the superior parietal pair to the
lowest value, which its own fibre counts contradict: 596 fibres of about 17
mm give the largest ratio, not the smallest. The mean fields of the same file
give a mean fibre count of 62.13 against the published 40.2 and reproduce
nothing.

The official mean-combination graph was also downloaded and aggregated with
the same confidence threshold. It gives mean adjacency `2.4398`, mean
weighted degree `66.4339`, and spectral gap `1.0608`, all far from the
published values (the retained reconstruction gives `1.5702`, `42.7547`, and
`0.7723`). A unilateral rather than bilateral entorhinal seed does not
recover the reported lobe separation. The audit is stored in
`results/connectome_weight_audit.csv`.

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
term distinguish the FEM model from the nodal reference. The reaction vector
is

```text
R_i(C) = int_Gamma alpha c_h (1 - c_h) phi_i ds,
```

integrated cell by cell with three-point Gauss-Legendre quadrature, exact
for its degree-four integrands, so no quadrature error enters the
comparison.

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
Panel (a) carries the anatomical inset of Fornari's figure 7: the connectome
with the four lobes in the lobe colours and every other vertex in grey. The
lobe membership mirrors `classify()` of the C++ test that computes the
curves, and the plot script asserts the same 58 four-lobe vertex count the
test enforces.

## Discretization study

The spatial study keeps `dt=0.4` and compares 1, 2, 4, and 8 P1 elements per
connection. The maximum difference between the four FEM biomarker curves and
the 8-element reference decreases as follows:

```text
elements per edge    DoFs    maximum difference [%]
1                     83       0.4992
2                   1213       0.1287
4                   3473       0.0246
8                   7993       0.0000 (reference)
```

The small negative transient observed with one element per edge disappears
from two elements per edge onward. This confirms that it is a coarse-mesh
effect and not a failure of the nonlinear solver.

The temporal study uses 8 elements per edge and compares the time steps
`dt=0.8, 0.4, 0.2, 0.1, 0.05`. Successive differences between biomarker curves give the rates

```text
dt       observed rate
0.4      1.0556
0.2      1.0262
0.1      1.0127
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
non-zero entries span 3.621 decades and 77.3% of them fall below 5% of the
maximum, so on a linear ramp three quarters of the connectome would render as
one indistinguishable colour. The two end values of the ramp are printed at
the bar. Absent connections are not on the scale at all; the 4629 zero
entries, the 83 diagonal ones among them, are drawn in one flat grey keyed on
the bar as "no connection", so a weak connection can never be confused with a
missing one.

Rows and columns keep the solver's own node order, with no permutation,
clustering or thresholding, and rows run from the bottom as in the printed
figure, so that the panel can be laid beside the published one and compared
cell for cell: the intra-hemisphere blocks sit in the lower-left and
upper-right quadrants, exactly where the paper describes them. Both colour
bars use the blue-to-red rainbow of the published figure. No cell is singled
out on the matrix: the panel's point is the heterogeneity of the network,
and the extreme pairs are verified against the paper in the domain
verification above instead.

The eight warm clusters along the diagonal, four per hemisphere, are the
connections internal to the cortical lobes, each bounded on the panel by a
thin dashed box computed from the region assignment. Every box is the
longest consecutive run of its lobe in the solver order, so it contains
cells of that lobe and of nothing else, and the script asserts as much
before drawing it: right frontal at indices 0-9, parietal 15-18, occipital
19-23, temporal 27-32, and the mirrored runs 41-50, 56-59, 60-64, 68-73 on
the left. Two vertices per hemisphere sit outside their lobe's box, because
the atlas enumeration interleaves the limbic belt, the cingulate between
the frontal and the parietal cortex and the entorhinal-parahippocampal pair
inside the temporal lobe: they are the paracentral gyrus, indices 10 and 51,
and the fusiform gyrus, indices 24 and 65, located by the strips. The mean
adjacency inside a lobe is between 1.91 and 4.25 times the mean of its
hemisphere block, which is what makes the clusters read warm on the
logarithmic ramp; these ratios are pinned by the verification harness. The
interleaving is also why the group strips break into 23 runs rather than
seven blocks: right hemisphere occupies indices 0 to 40, left 41 to 81, and
the brainstem 82.

For the public reconstruction used by the solver:

```text
unweighted degree range = 6 to 48
weighted degree range   = 2.0505 to 127.6435
adjacency range         = 0.0085 to 35.3221
```

The ranges coincide with the published 6-48, 2.1-127.6 and 0.01-35.32 at the
paper's printed precision, and the block structure shows the small-world
organization it describes. Per-region values are stored in
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
| mean adjacency | 1.57 | 1.5702 |
| mean weighted degree | 42.8 | 42.7547 |
| adjacency range | 0.01 to 35.32 | 0.0085 to 35.3221 |
| weighted degree range | 2.1 to 127.6 | 2.0505 to 127.6435 |
| model | `dc/dt = -L c + alpha c (1-c)`, equation (3.2) | same |
| conversion rate | `alpha = 0.5` | 0.5 |
| seed | entorhinal cortex at `c_0 = 0.1` | both entorhinal vertices |
| time integration | implicit, 100 steps of `dt = 0.4` years | same |
| biomarker | `C(t) = sum_I c_I(t)`, equation (4.1) | same, normalised to a percentage |

The paper defines the edge weight as `A_IJ = n_IJ / l_IJ`, the ratio of fibre
number and fibre length. Every aggregation of the public per-connection
fields was built and compared before one was retained:

```text
weighting of a fine connection, summed over parallel ones
                                      mean adjacency   mean w. degree
count_median / length_median               1.5702           42.7547   <- retained
electrical_connectivity_median             1.6669           45.3884
count_mean / length_mean                   2.2987           62.5910
sum(count_med) / mean(length_med)          1.3978           38.0599
paper reported values                      1.57             42.8
```

Only the first reproduces the published statistics, and it does so at printed
precision on every one of them, extremes and attaining regions included; the
`electrical_connectivity_median` field, the median of the per-subject ratios
rather than the ratio of the medians, misses by about six percent. The
match of eight independent statistics and four region identities rules out
coincidence, so the first weighting is the one used by the solver.

### The one quantity that is not reproduced

Figure 7 of the paper shows the four lobe biomarkers separated by roughly five
years, activating in the order temporal, frontal, parietal, occipital. Here
they cross the 50-percent level within `3.9e-7` years of each other.

Benchmark 23 resolves this. With `lambda_2 = 0.772` and `alpha = 0.5` the
Damkohler number `alpha / lambda_2` is `0.65`: the connectome homogenises
faster than the reaction grows, so coincident lobe curves are the correct
solution of equation (3.2) at unit weight scale. Reproducing the published
separation requires scaling the Laplacian by a factor near `0.005`. The
published figure therefore cannot come from the connectivity-weighted
Laplacian used at unit scale, whatever the reconstruction of the weights;
the exact reproduction of the published weights above makes this conclusion
sharper, since the reconstruction is no longer a candidate explanation.

The activation order is also worth recording. At any scaling that separates the
lobes, the frontal lobe is last here, because the frontal pole is the most
weakly connected node of the graph. The paper reports the same property in its
figure 9 discussion, where the frontal pole has the longest infection time of
all 83 seeding regions, so our ordering is consistent with the paper's own
network diagnostics even though it differs from its figure 7.

### Time discretization

The retained run uses `dt = 0.4` to match the paper. The temporal refinement
table in this benchmark shows the global 50-percent crossing moving from
`11.50` years at `dt = 0.8` to `13.71` years at `dt = 0.05`, so the headline
value carries a time-discretization error of order one year. The paper uses the
same step with a first-order implicit scheme and therefore carries a comparable
error; the agreement of the absolute crossing times should be read with that in
mind.

All CSV outputs in `results/` were regenerated from source in this working copy
and reproduce byte for byte.
