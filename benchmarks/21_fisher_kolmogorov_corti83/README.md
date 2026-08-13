# Benchmark 21: Deterministic Corti model on Budapest-83

## Goal

Apply the nonlinear metric-graph FEM solver to a deterministic version of the
Fisher-Kolmogorov model of Corti, Antonietti, Bonizzoni, Dede, and Quarteroni.
This is the final sequential application benchmark before performance,
reordering, and parallel-scaling studies.

Benchmark 17 is retained as the preliminary experiment on the 467-node
Budapest largest connected component. This benchmark instead uses the
83-region graph whose topology and weights were audited in benchmark 19.

## Literature correspondence

Corti et al. solve a patient-specific stochastic model on a 246-node
Brainnetome graph reconstructed from DWI and initialize it with patient PET
data. Those patient-specific graph and PET data are not available here.
Therefore this benchmark reproduces the published mathematical model, temporal
scheme, seven anatomical groups, and estimated mean reaction parameters, but
does not claim to reproduce the clinical trajectory of their patient.

The stochastic parameters are replaced by the reported estimated means:

| Anatomical group | alpha | Budapest-83 nodes |
|---|---:|---:|
| Frontal lobe | 0.1801 | 20 |
| Temporal lobe | 0.1421 | 14 |
| Parietal lobe | 0.0627 | 10 |
| Insular lobe | 0.1005 | 2 |
| Limbic lobe | 0.1351 | 16 |
| Occipital lobe | 0.0545 | 10 |
| Subcortical nuclei | 0.1147 | 11 |

The complete node-to-region assignment is stored in
`results/reaction_coefficients.csv`.

## Mathematical model

On each metric edge, the concentration satisfies

```text
partial_t c - partial_s(D_e partial_s c) = alpha(s) c (1-c).
```

Continuity and weighted Kirchhoff flux balance hold at every graph vertex. The
P1 finite element semi-discretization is

```text
M C_dot + K_D C = R(C).
```

The edge diffusivity is the normalized Budapest connectivity,
`D_e = w_e/max(w)`. The seven nodal reaction values are linearly interpolated
on every metric edge, and the reaction term is consistently integrated with
three-point Gauss-Legendre quadrature.

## Time discretization

The implementation follows the semi-implicit Crank-Nicolson construction used
by Corti et al. Diffusion and the first concentration factor are centered in
time, whereas the saturation factor is extrapolated with
`C* = 3/2 C^n - 1/2 C^(n-1)`, started from `C^(-1) = C^0`. The full update,
equation (4) of the paper on the metric graph, is

```text
(M + dt/2 K_D - dt/2 G(C*)) C^(n+1) = (M - dt/2 K_D + dt/2 G(C*)) C^n,
G(C*)_ij = int_Gamma alpha (1 - c*) phi_i phi_j.
```

Since `G(C*)` changes at every step, the matrix on the left is rebuilt and
factorized at every step; benchmark 22 measures this as the dominant cost.

```text
T                 = 20 years
dt                = 0.2 years
elements per edge = 4
DoFs              = 3473
```

The paper uses `dt=0.2` for calibration and `dt=0.02` for its stochastic
forward analysis. Here `dt=0.2` is appropriate for the deterministic
sequential application; temporal consistency of the scheme was already
verified in benchmark 17 and the nonlinear solver in benchmarks 18-20.

## Initial condition

Because the patient PET projection is unavailable, a documented synthetic
Alzheimer-like datum is used. Entorhinal and hippocampal regions start from
`c=0.10`; all remaining regions start from `c=0.01`. This choice creates
four seed nodes and must not be interpreted as patient data.

## Results

The solution remains inside the physical interval throughout the run. At
`T=20`, the full FEM field has range

```text
0.0238645 <= c_h <= 0.653239.
```

The mean concentration sampled at the 83 anatomical vertices increases from
`0.0143373` to `0.157989`. Every regional average grows monotonically.
The final regional means are:

```text
frontal       0.202253
temporal      0.167108
parietal      0.096674
insular       0.129048
limbic        0.210617
occipital     0.072868
subcortical   0.127734
```

The high frontal and low occipital values are consistent with their respective
reaction coefficients. The limbic curve is additionally elevated by the
synthetic seed. Corti et al. report high frontal and parietal patient values;
an exact ranking is not expected here because their PET initial condition,
Brainnetome topology, volumes, and calibrated patient graph are different.



## Visualization semantics

In Corti et al., marker size in the PET scatter plot is proportional to the
anatomical node volume, and marker size in the calibration comparison is
proportional to lobe volume. In Fornari et al., enlarged colored spheres are
also used to highlight selected seed regions; degree and edge color encode
connectivity in other network figures. Marker radius is therefore a
visualization attribute, not an additional numerical unknown.

Reliable regional volumes are not contained in the public Budapest files used
here. Our ParaView output consequently uses uniform node size and encodes
concentration by color. Assigning Corti-like variable radii without the
corresponding volume data would be misleading.

## Stored output

Report-ready results are stored in `results/`:

- `regional_averages.csv`: all seven regional curves and bounds;
- `reaction_coefficients.csv`: auditable anatomical assignment;
- `regional_averages.png` and `.pdf`: regional evolution.

The complete VTK time series remains generated under

```text
output/fisher_kolmogorov/corti83/solution.pvd
```

and contains `c`, `alpha`, and edge `diffusion`.

## Reproduce

From the project root:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target test_fisher_kolmogorov_corti83
./build-release/test_fisher_kolmogorov_corti83
python3 scripts/plot-fisher-kolmogorov-regions.py \
  output/fisher_kolmogorov/corti83/regional_averages.csv \
  --output output/fisher_kolmogorov/corti83/regional_averages.png
```

To visualize the graph in ParaView, open `solution.pvd` together with
`data/connectome/anatomy/brain_surface.vtk`. Display the surface in light
gray with opacity around `0.2`, color the graph by `c`, and keep a fixed
color range `[0,1]` across time.

## Literature verification

Checked against Corti, Antonietti, Bonizzoni, Dede and Quarteroni, sections 2.1
and 4, tables 1 and 3.

Confirmed against the paper:

- **Reaction coefficients.** All seven estimated means of table 1 are used
  verbatim: frontal `0.1801`, temporal `0.1421`, parietal `0.0627`, insular
  `0.1005`, limbic `0.1351`, occipital `0.0545`, subcortical `0.1147`.
- **Time discretization.** Equation (4) of the paper is implemented exactly:
  Crank-Nicolson on the diffusion and on the first concentration factor, with
  the saturation factor extrapolated as `3/2 c^k - 1/2 c^(k-1)`. The paper
  states `dt = 2e-1` for its simulations, which is the step used here; its
  `dt = 0.02` belongs to the stochastic forward analysis.
- **Horizon.** Twenty years, as in the paper's forward study.

Genuine deviations, all forced by data availability:

- The paper solves on a patient-specific 246-node Brainnetome graph from DWI;
  this benchmark uses the audited 83-region Budapest graph of benchmark 19.
- The paper initialises from a PET-PiB projection at patient age 61; this
  benchmark uses a documented synthetic Alzheimer-like seed. This is why the
  absolute concentrations here are far lower than the paper's.
- The paper does not state a normalisation of the edge weights. The choice
  `D_e = w_e/max(w)` is this project's, and benchmark 23 shows it is what
  places this benchmark in the regime where regional curves separate at all.

### Regional ranking against table 3

The absolute levels are not comparable, but the ordering of the seven regions
is. Table 3 of the paper reports the expected concentration per lobe at
`t = 20`; the final column below is this benchmark:

| Lobe | Corti et al., `t = 20` | rank | this benchmark | rank |
|---|---:|---:|---:|---:|
| frontal | 0.9289 | 1 | 0.2023 | 2 |
| limbic | 0.8905 | 2 | 0.2106 | 1 |
| temporal | 0.8699 | 3 | 0.1671 | 3 |
| insular | 0.8558 | 4 | 0.1290 | 4 |
| subcortical | 0.8413 | 5 | 0.1277 | 5 |
| parietal | 0.7738 | 6 | 0.0967 | 6 |
| occipital | 0.7336 | 7 | 0.0729 | 7 |

Twenty of the twenty-one pairwise orderings agree, giving a rank correlation of
`0.96`. The single inversion is frontal against limbic, and it is explained by
the initial condition: the synthetic seed places `c = 0.10` on the entorhinal
and hippocampal vertices, all of which are limbic, so the limbic curve starts
from `0.0325` instead of `0.01` and keeps that lead. Everything below the top
two is reproduced in the published order using nothing but the seven table-1
coefficients and the connectome topology.

Note that the prose of the paper describes "a higher concentration of misfolded
proteins inside the parietal and frontal lobes", whereas its own table 3 places
parietal second from last. The table is used here because it is unambiguous.
