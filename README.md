# FEMG

Lightweight C++ finite element code for parabolic PDEs on metric graphs.
The current implementation uses Boost Graph Library for graph topology and
Eigen for sparse linear algebra.

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Reproducible Benchmarks

Each report-ready benchmark can be regenerated independently. For example:

```bash
scripts/run-benchmark.sh 04_star_sine
```

To regenerate the full benchmark sequence:

```bash
scripts/run-all-benchmarks.sh
```

The scripts compile the required executable, run it, and copy the selected
outputs into the corresponding `benchmarks/*/results/` folder.

## Run the interval example

```bash
./build/test_heat
```

This writes ParaView files to:

```text
output/visualization/solution.pvd
output/visualization/solution_*.vtp
```

Open `output/visualization/solution.pvd` in ParaView.

## Four-pointed star example

The first parabolic test follows the four-pointed star geometry used in the
elliptic experiments of Abbate and Di Primio. The problem is

```text
u_t - Delta u = 0,    u(0) = 5,
```

so the exact solution is the constant function `u = 5`.

```bash
./build/test_heat_star_constant
```

This writes ParaView files to:

```text
output/visualization/star_constant/solution.pvd
output/visualization/star_constant/solution_*.vtp
```

Open `output/visualization/star_constant/solution.pvd` in ParaView.

The report-ready benchmark copy is stored in:

```text
benchmarks/02_star_constant/
```

## Four-pointed star with linear solution

The second parabolic star benchmark follows the linear-function test of Abbate
and Di Primio. The problem is

```text
u_t - Delta u = 0,    u(0,x,y) = x + y,
```

with Dirichlet values imposed at the four external vertices. The exact solution
is stationary:

```text
u(t,x,y) = x + y.
```

```bash
./build/test_heat_star_linear
```

This writes ParaView files to:

```text
output/visualization/star_linear/solution.pvd
output/visualization/star_linear/solution_*.vtp
```

The report-ready benchmark copy is stored in:

```text
benchmarks/03_star_linear/
```

## Four-pointed star with sine solution

The third parabolic star benchmark follows the sine-function test of Abbate and
Di Primio. The problem is

```text
u_t - Delta u = f,
```

with homogeneous Dirichlet conditions at all original vertices. The exact
stationary solution is

```text
u(t,x,y) = sin(2*pi*x),
```

and therefore

```text
f(x,y) = 4*pi^2*sin(2*pi*x).
```

```bash
./build/test_heat_star_sine
./build/test_heat_star_sine_convergence
```

Visualization output is written to:

```text
output/visualization/star_sine/solution.pvd
output/visualization/star_sine/solution_*.vtp
```

The report-ready benchmark copy is stored in:

```text
benchmarks/04_star_sine/
```

The spatial convergence table is stored in:

```text
benchmarks/04_star_sine/results/space_convergence.csv
```

## Four-pointed star with radial decay and reaction

This benchmark is a genuinely time-dependent parabolic test with a linear term
in the unknown:

```text
u_t - Delta u + r u = 0.
```

On each edge, with local coordinate `s` from the center to the external vertex,
the exact solution is

```text
u(t,s) = exp(-(pi^2/4 + r)t) cos(pi s / 2),
```

with `r = 1`. Homogeneous Dirichlet conditions are imposed at the four external
vertices, while the central vertex satisfies the natural Kirchhoff condition.

```bash
./build/test_heat_star_radial_decay
./build/test_heat_star_radial_decay_time_convergence
```

Visualization output is written to:

```text
output/visualization/star_radial_decay/solution.pvd
output/visualization/star_radial_decay/solution_*.vtp
```

The report-ready benchmark copy is stored in:

```text
benchmarks/05_star_radial_decay/
```

The temporal convergence plot and table are stored in:

```text
benchmarks/05_star_radial_decay/results/time_convergence.csv
benchmarks/05_star_radial_decay/results/time_convergence.pdf
```

## Eigenmode relaxation on the four-pointed star

This benchmark uses an exact solution built from an eigenmode of the graph
Laplacian:

```text
u_t - mu Delta u = 0.
```

The initial condition is asymmetric: one branch is hotter than the other
three. The nonconstant eigenmode decays in time and the solution relaxes toward
a constant state.

```bash
./build/test_heat_star_localized
```

Visualization output is written to:

```text
output/visualization/star_localized/solution.pvd
output/visualization/star_localized/solution_*.vtp
```

The report-ready benchmark copy is stored in:

```text
benchmarks/06_star_localized/
```

## Fisher-Kolmogorov validation

The nonlinear solver is first validated on a one-dimensional sensitivity test
with 200 P1 elements, following the setup of Weickenmeier et al. Both the
semi-implicit method of Corti et al. and fully implicit Backward Euler with
Newton are available.

```bash
./build/test_fisher_kolmogorov_logistic
./build/test_fisher_kolmogorov_1d_sensitivity
```

The report-ready profiles, figure, parameters, and reproduction commands are
stored in:

```text
benchmarks/18_fisher_kolmogorov_1d_sensitivity/
```

The verified 83-region Budapest/FreeSurfer application compares the nodal
network model from the literature with the P1 metric-graph FEM formulation:

```bash
python3 scripts/prepare-fornari-connectome.py
./build/test_fisher_kolmogorov_fornari83
```

Its complete documentation and report-ready biomarker comparison are stored in:

```text
benchmarks/19_fisher_kolmogorov_fornari83/
```

The conversion-rate sensitivity study reproduces the deterministic
Fisher-Kolmogorov parameter trend for `alpha=0,0.1,...,0.5` and checks the
conserved nodal average and FEM metric integral:

```bash
python3 scripts/audit-fornari-connectome.py
python3 scripts/plot-fisher-alpha-sensitivity.py
```

Report-ready output is stored in:

```text
benchmarks/20_fisher_kolmogorov_alpha_sensitivity/
```

## Graph input format

The reader accepts the original edge-list format:

```text
n_vertices n_edges
u v length n_cells
...
```

It also accepts a coordinate-aware format:

```text
n_vertices n_edges
x0 y0
x1 y1
...
u v length n_cells
...
```

The coordinate-aware format is used for planar visualization in ParaView.
Internally, global DoFs are numbered as in the extended graph construction:
artificial edge nodes first, edge by edge, and original graph vertices last.

## Validation

Spatial convergence:

```bash
./build/test_heat_convergence
python3 scripts/plot-convergence.py output/convergence/space_convergence.csv
```

Temporal convergence:

```bash
./build/test_heat_time_convergence
python3 scripts/plot-time-convergence.py output/convergence/time_convergence.csv
```

Generated convergence files are written to:

```text
output/convergence/
```

## Graph setup diagnostic

```bash
./build/test_graph_setup
```

This prints the global DoF numbering associated with the graph input file.


The deterministic heterogeneous extension using the seven regional reaction
means reported by Corti et al. is documented in:

```text
benchmarks/21_fisher_kolmogorov_corti83/
```


The optimized single-thread baseline preceding graph reordering and
parallelization is documented in:

```text
benchmarks/22_fisher_kolmogorov_sequential_performance/
```

The diffusion-scaling study that explains why the lobe biomarkers of benchmark
19 coincide while those of benchmark 21 separate, and that establishes the
validity boundary of the consistent-mass P1 discretization, is documented in:

```text
benchmarks/23_fisher_kolmogorov_diffusion_scaling/
```

```bash
python3 scripts/study-fisher-diffusion-scaling.py \
  --output-dir output/fisher_kolmogorov/diffusion_scaling
python3 scripts/plot-fisher-diffusion-scaling.py \
  output/fisher_kolmogorov/diffusion_scaling \
  --output benchmarks/23_fisher_kolmogorov_diffusion_scaling/results/diffusion_scaling.png
```

`test_fisher_kolmogorov_fornari83` accepts, after the output directory, the
conversion rate, the final time, a uniform scaling of the connectivity weights,
and the time scheme (`be`, `cn` or `nodal`):

```bash
./build/test_fisher_kolmogorov_fornari83 1 0.4 output/run 0.5 40 1.0 be
```

The connectome domain, the regional seeding vulnerability and the anatomical
progression of the spreading are documented in:

```text
benchmarks/24_connectome_topology/
benchmarks/25_connectome_seeding_vulnerability/
benchmarks/26_connectome_progression/
```

```bash
python3 scripts/plot-connectome-regions.py \
  --output-dir benchmarks/24_connectome_topology/results
python3 scripts/study-fisher-seeding-vulnerability.py \
  --output-dir output/fisher_kolmogorov/seeding
python3 scripts/plot-connectome-progression.py \
  --output benchmarks/26_connectome_progression/results/anatomical_progression.png
```

The figures of these three benchmarks share `scripts/connectome_style.py`,
which reads the anatomical group of every vertex from the CSV written by
`test_fisher_kolmogorov_corti83` so that the figures and the solver cannot
disagree about the classification.

## Report figures

To stage every report-ready figure of the connectome chapter:

```bash
scripts/collect-report-figures.sh
```

This copies the PDFs into `report/images/`, which is not tracked because the
figures are already stored with their benchmarks.

## Plotting scripts

The plotting scripts need `matplotlib` and `numpy`. The Python shipped with the
`gcc-glibc` toolchain has an incomplete `matplotlib` installation, so run them
with the system interpreter:

```bash
/usr/bin/python3 scripts/plot-fisher-fornari83.py --help
```
