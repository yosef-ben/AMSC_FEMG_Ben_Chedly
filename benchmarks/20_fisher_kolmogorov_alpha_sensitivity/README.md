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
every positive `alpha` reaches `99.99%` or more of the stable state within the
horizon, and the 50-percent crossing time grows from `11.1` years at
`alpha = 0.5` to `59.3` years at `alpha = 0.1` in the nodal reference. The
`alpha = 0` case conserves the nodal average to `3.9e-16` and the FEM metric
integral to `9.7e-17`.

The nodal and FEM curves differ in amplitude because a value prescribed at a
graph vertex has a connectivity-dependent metric support in the FEM domain.
This is a model difference, not a discretization error, and it is left
uncorrected.
