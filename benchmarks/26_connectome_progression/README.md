# Benchmark 26: Anatomical Progression

## Goal

Show the spreading itself, on the anatomy, rather than through an integrated
biomarker, in the visual grammar the reference works use for that purpose.

## The two figures

`results/anatomical_progression` is a **mid-sagittal section** at three stages.
`results/activation_time` is a single image carrying the whole history: every
degree of freedom coloured by the year its concentration reaches `0.5`.

Both come from the deterministic Corti model of benchmark 21, refined to eight
elements per connection, over a horizon of 60 years, and are read from the
`solution_*.vtp` files written by `test_fisher_kolmogorov_corti83`.

## How the section is built

Weickenmeier et al. show an opaque brain cut at the mid-sagittal plane with
their tetrahedral mesh visible in neutral grey where the concentration is low.
Three things follow from copying that construction rather than its palette:

- The pial surface is clipped to the half-space `x >= 0` and the graph to
  `x <= 0`, on the same plane. The surface is watertight, `vtkFeatureEdges`
  reports zero boundary and zero non-manifold edges, so the cut is a genuine
  section and its cap is drawn slightly darker to read as a cut face. One depth
  buffer then occludes correctly: the network sits inside the anatomy instead
  of floating over it, and no layer trick is used.
- The zero of the colour ramp **is** the neutral grey of the discretisation, so
  the network reads as a pale wireframe where the concentration is low and the
  lesion emerges from it. This is what the previous versions of this figure got
  wrong: `inferno` puts black at `c = 0`, so the early stages were a dark mass.
  Nothing is thresholded away; the whole field is drawn continuously.
- Tube and vertex radii are constant, `0.55` mm and `2.6` mm. Thickness encodes
  nothing; only colour carries a value.

White marks are the exact `c = 0.5` level set, computed segment by segment on
the P1 interpolant: on a segment with endpoint values `a` and `b` straddling
the level, the crossing sits at `(0.5 - a)/(b - a)` along it. Their scatter is
itself a result, discussed below.

## How the stages are chosen

By a rule, not by eye: the three stored times at which the network mean
concentration reaches `0.2`, `0.4` and `0.6`. That is the same quantity the
reference plots as its biomarker, so the staging is a property of the solution.

```text
stage   t [yr]   mean c   network length above c = 0.5   level-set marks in section
1        23      0.2069                            1.9%                          40
2        31      0.4063                           22.8%                         223
3        38      0.6018                           72.3%                         247
```

The length fraction is exact for a P1 field: on a segment with endpoint values
`a`, `b` the supra-threshold fraction is `0`, `1`, or `(max - 0.5)/(max - min)`,
weighted by segment length.

Equally spaced times do not work and it is worth recording why. At `t = 15`
years no degree of freedom anywhere exceeds `0.454`, so the first panel is a
blank wash; by `t = 45` years more than 95% of the network length is above
`0.5`, so the last panel is a uniform mass.

## What the activation-time figure settles

The obvious reading of a spreading figure is that a front travels outward from
the seed. On this connectome that reading is wrong, and the second figure shows
why. Over the 83 region vertices, the activation time correlates

```text
with the local reaction coefficient   r = -0.841
with the graph distance from the seed r = +0.072
```

and the graph eccentricity from the seed is `2`: every region lies within two
connections of the entorhinal and hippocampal seed. The order in which the
network is recruited is therefore set by the regional conversion rates,
`0.1801` in the frontal group against `0.0545` in the occipital one, and not by
any travelling front. This is the same conclusion benchmark 23 reaches from the
Damkohler number, arrived at independently.

Two consequences. The anterior-to-posterior gradient visible in both figures is
a **reaction-rate pattern**, not a propagation pattern, and the caption says
so. And no arrows are drawn: an arrow would assert a direction of travel that
the solution does not contain. The scattered `c = 0.5` marks in the section are
the visual counterpart of the same fact, the level set is not a localised
surface.

`10` of the `10170` degrees of freedom never reach `0.5` within the horizon;
they take the late end of the ramp and are counted in the caption.

## Two deliberate choices

The horizon is 60 years, not the 20 years of benchmark 21. Corti et al. stop at
20 because that is their clinical question; with their mean coefficients the
concentration only reaches `0.61` by then. The model, the coefficients, the
seed and the time step are unchanged, and the run is stored separately, so
benchmark 21 is untouched.

Only the half-space `x <= 0` of the network is drawn. This is justified by
bilateral symmetry, the mean absolute left-right difference over the 41 mirror
region pairs stays below `0.011`, under 4% of the mean concentration at every
stage, but it is still a sub-picture of the run.

## What cannot be reproduced, and is not faked

The reference fills its brain silhouette because its unknown lives on a
tetrahedral volume mesh. Ours lives on 1130 one-dimensional connections, so the
coloured structure occupies a lens-shaped core and there is bare tissue around
it. That rim is an absence of degrees of freedom, not an absence of disease.
Filling it, by interpolating the graph field into the volume or by assigning
pial triangles to their nearest region, is the one move that would make the
figure look like theirs, and it is exactly the fabrication these figures avoid.

Two further points a reader will otherwise misread. Pale grey means low, never
healthy: the initial condition is `c = 0.01` everywhere and the network minimum
rises monotonically, so no part is ever at zero, whereas the reference's white
tissue genuinely is. And the interior nodes of a connection lie on the straight
chord between two region centroids; that is the discretization the solver used,
not an anatomical fibre path.

## Reproduce

Run the commands in `commands.txt` from the project root.
