# Benchmark 24: Topology of the 83-Region Connectome

## Goal

Document the domain used by every connectome experiment in this project, in
the form used by the reference works, so that the graph can be inspected
independently of any simulation.

## Content

`results/connectome_regions` shows the 83 vertices at their anatomical
coordinates, coloured by the seven anatomical groups of Corti et al., first
without and then with the 1130 connections, in a sagittal and an axial view.
Line width is proportional to the connectivity weight.

`results/connectome_connectogram` shows the connections above 5% of the
strongest one on a circular layout grouped by region, which is the threshold
Corti et al. use for their connectogram, together with the connectivity summed
over each pair of groups.

## Relation to the reference figures

Corti et al. present the same object in four panels: the MRI brain surface, the
DWI tractography, the graph vertices and the graph. The first two are not
reproduced. The tractography is not available, the graph here is taken from the
public Budapest Reference Connectome rather than reconstructed from patient
images, and drawing an MRI surface behind our vertices would suggest a
provenance the data does not have.

The region assignment is read from `reaction_coefficients.csv`, the file
written by `test_fisher_kolmogorov_corti83`, so the figure and the solver can
never disagree about which vertex belongs to which group.

## What the connectivity matrix shows

The aggregated matrix is dominated by short-range connections: frontal to
parietal reaches `183` and parietal to occipital `202`, while frontal to
temporal is only `1`, spread over 14 connections whose individual weights do
not exceed `0.19`. This is a direct consequence of the weight definition, the
fibre count divided by the fibre length, which penalises exactly the long
association bundles that join the frontal and temporal lobes. It is a property
of the reconstruction, not of the solver, and it is worth keeping in mind when
reading the regional results of benchmarks 21 and 25.

## Reproduce

```bash
python3 scripts/plot-connectome-regions.py \
  --output-dir benchmarks/24_connectome_topology/results
```
