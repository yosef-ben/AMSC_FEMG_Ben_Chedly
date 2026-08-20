# Benchmark 24: Topology of the 83-Region Connectome

## Goal

Document the domain used by every connectome experiment in this project, in
the form used by the reference works, so that the graph can be inspected
independently of any simulation.

## Content

`results/connectome_regions` shows the graph inside the pial surface in the
sagittal, coronal and axial projections. The top row colours the 83 vertices by
the seven anatomical groups of Corti et al., with the colours of the legend
of their figure 2, sampled from the published panels, so the two figures can
be read against each other group by group. The choice trades colour-vision
robustness for comparability: under severity-one dichromacy simulation the
smallest pairwise OKLab (x100) separation is 3.3, insular against limbic in
tritanopia, against 16.5 in normal vision; the renders mitigate this with
dark marker edges and lightness differences. The same palette is used
wherever the seven groups appear: here, in the regional curves of benchmark
21, in the seeding study of benchmark 25 and in the group strips of the
adjacency panel of benchmark 19.
The bottom row colours and sizes the 1130 connections by their connectivity
weight on the rainbow ramp of the brain network figure of Fornari et al.,
with the uniform silver spheres of that figure at the vertices. Colour and
width follow the square root of the weight, a mapping stated here because
the references declare none: the distribution is strongly skewed, a linear
ramp paints nearly everything at the cold end, and the square root keeps the
mid-range readable while the short association bundles keep the warm
colours. The logarithmic reading of the same weights lives in the adjacency
panel of benchmark 19, where cell-level comparison is the point.

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
over each pair of groups. The arcs follow the colouring of their figure 3:
each connection is a gradient between the colours of the two groups it joins,
thickness grows with the weight, and the stronger connections are drawn on
top. The group matrix keeps a perceptual ramp with its values printed in the
cells, leaving the rainbow to the renders that reproduce reference figures.

`results/connectome_views` reproduces the brain network figure of Fornari et
al., their figure 5, in its own composition: the two-by-two arrangement of a
sagittal, a coronal, a longitudinal and an unlabelled oblique view, the last
showing the graph without the surface, as the printed panel does; uniform
silver spheres at the 83 vertices, and the 1130 connections drawn with colour
and thickness both growing with the fibre number, over a rainbow bar
labelled only at its ends, centred under the views. Two departures from the printed figure are
deliberate and stated in the report caption as well. The vertex coordinates
and the pial surface are the ones distributed with the public Budapest
Reference Connectome viewer, not the MRI-based brain model of the reference,
so the layout matches anatomy but not their rendering point for point. And
the reference does not declare how fibre number maps to colour and width, so
the mapping here is a stated choice: both follow the square root of the
fibre number, between 1 and 595.5. The distribution is strongly skewed, so
the strongest short association bundles, the superior parietal to precuneus
pair above all, carry the warm colours while the square root keeps the
mid-range distinguishable; the printed figure shows the same broad spread of
warm edges, which a purely linear ramp of so skewed a distribution could not
produce. The
rainbow ramp breaks the colour conventions used everywhere else in this
project; it is kept in these network renders, here and in the weight row of
`connectome_regions`, because matching the reference's own colouring is what
makes ours and theirs directly comparable.

`results/lobe_connectivity` sums the same 1130 connections over the four
cortical lobes of Fornari et al., the partition of every biomarker curve of
the report, plus the 25 remaining regions (insular, limbic, subcortical and
the brainstem): one line per pair of groups with the total on it, the width
growing with the square root of the value, the entorhinal seed starred. It
is the quotient graph behind the activation order of the report: the
temporal lobe leads to the parietal (75.0) and occipital (55.6) lobes, while
the frontal lobe hangs on the parietal lobe (190.5) and the deep regions
(220.3), its direct couplings to the temporal and occipital lobes being 1.0
and 2.3. These sums use the four-lobe partition and cannot be compared cell
by cell with the seven-group matrix of `connectome_connectogram`.

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
