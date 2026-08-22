# Benchmark 20: Conversion-Rate Sensitivity on the Connectome

## Goal

Study the effect of the Fisher-Kolmogorov conversion rate on the 83-region
connectome. This is the direct deterministic analogue of the conversion-rate
experiment in Fornari et al. and tests the nonlinear application beyond a
single parameter choice.

## Problem

The nodal reference solves

```text
c_dot + L c = alpha c (1-c).
```

The metric-graph formulation solves

```text
M C_dot + H C = R_alpha(C).
```

The two entorhinal regions are initialized at `c=0.1` and all other
regions at zero. We use Backward Euler with `dt=0.4` up to `T=80` and
compare

```text
alpha = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5.
```

One P1 element per connection is used here because the purpose is a direct
comparison between the 83-unknown nodal model and its lowest-order
metric-graph counterpart. Spatial and temporal refinement were verified
separately in benchmark 19.

## Expected behaviour

Increasing `alpha` must move the sigmoid biomarker curve to earlier times.
For `alpha=0` there is diffusion without conversion. The nodal model
conserves the average over the 83 regions, while the FEM model conserves the
normalized metric integral `1^T M C / (1^T M 1)`. Every positive conversion
rate eventually drives the solution toward the stable state `c=1`.

## Reproduce

Run the commands in `commands.txt` from the project root.

## Stored results

`results/alpha_sensitivity.csv` reports threshold-crossing times, final
global concentrations, concentration bounds, and both conservation checks.
The nodal and FEM curves start from different global averages because a value
prescribed at a graph vertex has a connectivity-dependent metric support in
the FEM domain. No amplitude rescaling is applied to hide this model
difference.
`results/alpha_sensitivity.png` and the PDF version are ready for the
report.

## Use in the report

The figure is a record and is not in the report, which keeps one sentence
of this study where the regional rates are introduced: lowering the rate
from 0.5 to 0.1 moves the 50-percent crossing of the finite element network
from 12.68 to 67.81 years while every positive rate still drives the network
to the misfolded state.

## Interpretation

This benchmark reproduces the qualitative conclusion of the reference
conversion-rate study: reducing `alpha` delays the transition but does not
create an intermediate stable state when `alpha>0`. It is not a fit to
clinical data, and the nodal/FEM difference is retained as a numerical-model
comparison.

## Literature verification

Checked against Fornari, Schafer, Jucker, Goriely and Kuhl, *Prion-like
spreading of Alzheimer's disease within the brain's connectome*, J. R. Soc.
Interface 16 (2019) 20190356, section 5.1 and figures 10 and 11.

The experiment is the deterministic analogue of the paper's conversion-rate
study, and the parameter set matches it exactly: `alpha = 0.5, 0.4, 0.3, 0.2,
0.1, 0` over a horizon of 80 years, seeded in the entorhinal cortex at
`c_0 = 0.1`.

The paper's conclusion is that "irrespective of the conversion rate alpha, the
misfolded protein concentration of the Fisher-Kolmogorov model always converges
towards the fully misfolded state with a biomarker abnormality of C = 100%",
and that lowering `alpha` only delays the transition. Both are reproduced:
every `alpha >= 0.2` reaches `99.98%` or more of the stable state within the
horizon, `alpha = 0.1` is still rising at 80 years (`88.6%` nodal, `77.1%`
FEM, the metric average) and the 50-percent crossing time grows from `11.1`
years at `alpha = 0.5` to `59.3` years at `alpha = 0.1` in the nodal
reference, from `12.7` to `67.8` years in the FEM. The product of the
conversion rate and the FEM crossing time stays between 6.3 and 6.8 years
over the five positive rates, the scaling of a logistic growth from the
seed level. The `alpha = 0` case conserves the nodal average to `1.1e-14`
and the FEM metric integral to `8.0e-16` (`alpha_sensitivity.csv`).

The FEM runs use the consistent mass at `rho = 1`, where benchmark 23 shows
it valid: the Damkohler numbers of the five positive rates run from 1.6 to
7.9 (`Da_lobe` from 0.09 to 0.45 with the lobe rate 1.107 of benchmark 23),
the transient minimum stays above `-1.4e-4` and every lobe is synchronized.
The FEM biomarker plotted is the metric average `1^T M C / (1^T M 1)`, the
quantity the scheme conserves at `alpha = 0`; its 50-percent crossing at
`alpha = 0.5`, 12.679 years, agrees with the vertex-average crossing of
benchmark 19, 12.682 years, to three thousandths of a year.

The nodal and FEM curves differ in amplitude because a value prescribed at a
graph vertex has a connectivity-dependent metric support in the FEM domain.
This is a model difference, not a discretization error, and it is left
uncorrected.
