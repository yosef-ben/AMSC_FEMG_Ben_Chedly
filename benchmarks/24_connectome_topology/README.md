# Benchmark 24: Topology of the 83-Region Connectome

## Goal

Document the domain used by every connectome experiment in this project, in
the form used by the reference works, so that the graph can be inspected
independently of any simulation.

## Content

`results/connectome_regions` shows the graph inside the pial surface in the
sagittal, coronal and axial projections. The top row colours the 83 vertices by
the seven anatomical groups of Corti et al.; the bottom row colours the 1130
connections by their connectivity weight, as in the brain network figure of
Fornari et al.

The panels are rendered with VTK, the engine ParaView is built on, so the
figures and an interactive ParaView session show the same picture. The surface
is the same `brain_surface.vtk` opened in ParaView, in the same coordinate
frame as the vertices: nothing is rescaled or repositioned to make the two
agree, and all 83 vertices fall inside it. Regions are drawn as shaded spheres
and connections as shaded tubes, with depth peeling for the translucent
surface. Every panel uses one camera scale, so the three projections are
directly comparable.

`results/connectome_connectogram` shows the connections above 5% of the
strongest one on a circular layout grouped by region, which is the threshold
Corti et al. use for their connectogram, together with the connectivity summed
over each pair of groups.

`results/connectome_views` reproduces the brain network figure of Fornari et
al., their figure 5, in its own composition: the two-by-two arrangement of a
sagittal, a coronal, a longitudinal and an unlabelled oblique view, uniform
silver spheres at the 83 vertices, and the 1130 connections drawn with colour
and thickness both growing with the fibre number, under a rainbow bar
labelled only at its ends. Two departures from the printed figure are
deliberate and stated in the report caption as well. The vertex coordinates
and the pial surface are the ones distributed with the public Budapest
Reference Connectome viewer, not the MRI-based brain model of the reference,
so the layout matches anatomy but not their rendering point for point. And
the reference does not declare how fibre number maps to colour and width, so
the mapping here is a stated choice: both scale linearly between the smallest
and the largest fibre number, 1 and 595.5. The distribution is strongly
skewed, so most connections stay thin and blue and a handful of short
association bundles, the superior parietal to precuneus pair above all,
carry the warm colours; the printed figure shows the same aspect. The
rainbow ramp breaks the colour conventions used everywhere else in this
project; it is kept in this one figure because matching the reference's own
colouring is what makes the two directly comparable.

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
parietal reaches `175` and parietal to occipital `190`, while frontal to
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
python3 scripts/plot-connectome-views.py \
  --output-dir benchmarks/24_connectome_topology/results
```
