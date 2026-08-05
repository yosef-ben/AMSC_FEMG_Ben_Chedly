# Benchmark 26: Anatomical Progression

## Goal

Show the spreading itself, on the anatomy, rather than through an integrated
biomarker. This is the figure the reference works use to argue that the model
reproduces the clinical staging of the disease.

## Content

Rows are simulation settings, columns are four instants of the same run. The
vertices carry the concentration at the anatomical coordinates of the region;
the connections carry the mean of their endpoints, which is exactly what the P1
field with one element per connection is. The colour range is fixed to `[0,1]`
in every panel, so panels can be compared directly.

The nodal reference is plotted throughout. It is the model for which this
figure exists in the reference works, and, unlike the consistent-mass P1
formulation, it remains bounded at every connectivity scaling shown here
(benchmark 23).

## Relation to the reference figures

Fornari et al. show the effect of lowering the conversion rate on the
anatomy. The first two rows are that experiment: at `alpha = 0.3` the same
concentration is reached roughly six years later than at `alpha = 0.5`, with an
unchanged spatial pattern.

Weickenmeier et al. show whole-brain progression in sagittal, transverse and
coronal sections. Those figures come from a continuum model on a tetrahedral
brain mesh with about 80000 unknowns. Our domain is a metric graph with 83
vertices, so a sectioned volumetric field cannot be produced from it and none
is drawn. The projection used here is a plain orthographic view of the graph;
no interpolation into the brain volume is performed.

## What the figure shows about the model

The third row is the reason the figure is worth including. At the connectivity
scaling of the paper the concentration is spatially uniform at every instant:
all vertices share the same colour, and the panels differ only in overall
level. A visible front requires the reduced scaling, and the third row shows
it: at `t = 6` and `t = 10` years the medial temporal seed region is clearly
ahead of the rest of the graph, and at `t = 18` the frontal pole, the weakest
connected vertex, is still behind.

Reporting the uniform rows alongside the structured one is deliberate. Showing
only the reduced-scaling row would suggest that the model of the reference
paper produces an anatomical front, which on this graph it does not.

## Reproduce

The figure reads the profiles produced by benchmarks 20 and 23, so run those
first. Then run the command in `commands.txt`.
