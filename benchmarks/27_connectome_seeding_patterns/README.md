# Benchmark 27: Seeding Patterns and Expected Progression

## Goal

Show, in the construction of figure 1, panel (c), of Fornari et al., what the
network model produces for the two clinical seedings of Alzheimer's disease,
and compare it with the progression the literature expects: tau inclusions
first in the transentorhinal region, spreading through the temporal lobe to
the interconnected neocortex (Braak staging); amyloid-beta deposits first in
the neocortex, descending towards the deep structures, with the subcortical
nuclei last (Weickenmeier et al., section 2.7, after Jucker and Walker).

## Setup

Both rows solve the metric-graph FEM with the lumped mass and the fully
implicit scheme (`be_lumped`), at one element per connection, with
`alpha = 0.5`, `dt = 0.4` years and `rho = 0.005`, the scaling at which
benchmark 23 obtains a separation of the magnitude published by Fornari et
al. (`Da` of order one hundred). The lumped mass is used because at this
scaling the consistent-mass FEM leaves the physical range and diverges, the
validity boundary documented in benchmark 23; the lumped solution stays in
`[0,1]` throughout and the plot script asserts it. The earlier nodal-model
version of this figure gave the same qualitative picture with stages at
10.4, 14.8, 18.8 and 0.8, 4.4, 8.8 years.

Seeding, at `c0 = 0.1`:

- tau: the two entorhinal cortices, the default of the executable;
- amyloid-beta: the 54 neocortical vertices of the four cortical lobes,
  the `neocortex` seeding mode of `test_fisher_kolmogorov_fornari83`; the
  entorhinal and parahippocampal cortices, temporal in the partition but
  allocortex, which amyloid-beta reaches only at its second stage (Thal
  phase 2), are not seeded.

The three stages are selected by a rule, not by eye: the first instants at
which the network mean reaches 10, 40 and 80 percent.

The scaling is the whole point and a control run shows it. `rho = 1`, the
transport scale of the references, gives `Da_lobe = 0.58` for this model and
the transport synchronizes the regions: the tau run repeated at `rho = 1`
(`results/tau_rho1_*.csv`) reaches its first stage at 8.8 years with every
vertex between 0.101 and 0.131, its lobes crossing 50 percent within 0.027
years of each other (12.705 to 12.732), so its three renders would be three
uniform mantles. The progression is visible only at `Da_lobe` well above
one. This is a calibration of the transport time scale against the
reference, not a biological input: with a uniform conversion rate, the only
biology the runs contain is the location of the seeds, and the figure tests
whether the connectivity orders the regions as the clinic does.

## Results

```text
tau      stages at 15.2, 21.6, 27.6 years
amyloid  stages at  1.2,  5.2, 12.0 years

amyloid-beta, mean group activation (first stored crossing of c = 0.5):
temporal 4.4, parietal 4.4, occipital 4.4, frontal 4.4,
insular 9.6, limbic 11.68, subcortical 13.05 years;
brainstem vertex 17.2 years, the last of the 83

tau, lobe 50-percent crossings: temporal 16.50, occipital 21.65,
parietal 22.95, frontal 25.87 years
```

The amyloid ordering is the expected one, the four amyloid phases of Thal in
sequence: the seeded neocortical lobes activate first and together, then
the insula and the limbic belt (the allocortex, entorhinal 11.6 and 12.4
years), then the subcortical nuclei, more than eight years behind the cortex,
and the brainstem last of all. With the earlier seed, which included the
entorhinal and parahippocampal cortices among the 58 vertices of the four
lobes, the stages were 0.8, 4.8 and 11.2 years and the limbic mean 9.78
years, pulled forward by its seeded members. For tau
the first activated vertices are the two entorhinal seeds and the spreading
reaches the temporal neighbourhood before the rest of the cortex, as in the
Braak sequence; the ordering of the late lobes differs from the clinical one
for the reason recorded in benchmark 23, the weak connectivity of the
frontal pole.

The occipital-before-parietal order is a property of the group means, not of
the direct couplings: the temporal lobe is coupled to the parietal one more
strongly than to the occipital one (75.0 against 55.6 in total connectivity)
and the first parietal region does cross first (19.4 against 19.6 years),
but the twelve parietal regions carry a slow tail (postcentral and
paracentral, 23.2 to 26.0 years) while the eight occipital ones lie together
(19.6 to 23.4 years), so the parietal mean crosses last. The numbers come
from `tau_profiles.csv` and the edge list; the verification recomputes them
and `results/lobe_crossings` draws them, one dot per region on one row per
lobe with the crossing of the lobe mean marked.

## The same seedings with the regional reaction rates

The uniform conversion rate is the only reason the late lobes activate in
the order of the connectivity. Repeating both seedings with the seven
regional rates of Corti et al. (`reaction_coefficients.csv` of benchmark 21,
read by the executable and rescaled to the same vertex mean 0.5, so the
reaction budget and both Damkohler numbers are unchanged and only the
distribution of the rate differs) recovers the clinical sequence:

```text
tau, lobe 50-percent crossings [years]
  uniform rate    temporal 16.50, occipital 21.65, parietal 22.95, frontal 25.87
  regional rates  temporal 14.70, frontal 19.91, parietal 29.71, occipital 32.39
  Fornari et al.  temporal, frontal, parietal, occipital

tau, lobe means at the three stages [percent]
  uniform   t=15.2  38.2 / 0.4 / 2.3 / 3.5      (temporal/frontal/parietal/occipital)
            t=21.6  86.1 / 12.6 / 35.8 / 49.4
            t=27.6  98.6 / 66.7 / 88.2 / 93.9
  regional  t=14.0  42.8 / 1.3 / 0.4 / 0.6
            t=20.4  90.6 / 57.2 / 6.3 / 5.7
            t=28.8  99.5 / 99.1 / 44.3 / 31.6

anterior-posterior centroid of the regions above 0.5, at the three stages
  uniform    +8.5 -> -14.6 ->  -7.2 mm     (the spreading moves posteriorly)
  regional   +2.2 ->  +8.9 -> +11.9 mm     (it moves anteriorly, as the
                                            clinical arrows do)

amyloid-beta, regional rates, group means of the first crossing [years]
  seeded neocortex 5.6, limbic 10.9, insular 11.2, remaining subcortical
  13.3, brainstem 18.4: the downward progression is preserved
  stages at 1.2, 5.6 and 13.2 years
```

The rates run from 0.2177 (occipital) to 0.7195 (frontal) after rescaling, a
factor 3.3, and the solution stays in `[0,1]`.

The gain is tau's alone, and the amyloid row is measurably worse. With the
uniform rate the seeded mantle rises as a block (spread between its 54
vertices 0.00, 0.02 and 0.01 at the three stages) so the only visible
contrast is the one to the deep structures, which is the Thal picture. With
the regional rates the mantle acquires an anterior-to-posterior gradient
that the staging drawings do not show (spread 0.10, 0.58 and 0.32; the
concentration-weighted centroid moves from -4.9 to -2.3 mm with the uniform
rate and from +3.9 to +9.8 mm at the middle stage with the regional one),
and the cortex-to-deep contrast at the last stage falls from +0.56 to +0.43.
The ordering survives, the seeded neocortex still crossing before the
allocortex, the deep nuclei and the brainstem, although the margin between
the last cortical group and the first allocortical one narrows from 5.2 to
1.0 years (the four seeded lobes cross at 3.2, 4.0, 8.4 and 9.9 years
instead of together at 4.4). The gradient cannot be judged against the
staging drawings, which shade the whole neocortex at the first phase and
resolve no ordering inside it, and it is the kind of detail these
coefficients carry: Corti et al. estimate them from amyloid-beta PET, so
applying them to the tau seeding, where they restore the clinical lobe
order, is a transfer across proteins. The report shows the tau row alone,
the amyloid study having been dropped from it with the motivation stated in
the staging section, and `results/seeding_patterns_regional_full_expected`
keeps both rows as the record behind these numbers. Below the render rows
the report places `results/regional_biomarker_curves`, the same regional tau
run read through the four lobe biomarker curves in the presentation of
figure 7 of Fornari et al., with the network mean, the 50-percent level and
the three stage instants marked.

## The control that removes the protein transfer

Because those coefficients are estimated for amyloid-beta, the tau result
must not depend on their particular values, and it does not. Repeating the
tau run with `results/synthetic_rate_field.csv`, a field linear in the
anterior-posterior coordinate of the vertices that shares with them only the
ratio 3.3 between its fastest and its slowest value (and, like them, is
rescaled by the executable to the vertex mean 0.5), gives

```text
temporal 17.46, frontal 21.12, parietal 26.75, occipital 31.25 years
```

the same sequence the coefficients of the reference produce (14.70, 19.91,
29.71, 32.39) and the one the uniform rate cannot produce (16.50, 21.65,
22.95, 25.87).

`scripts/study-rate-permutations.py` then asks what a regional field must
satisfy, by assigning the seven values to the seven groups in all 5040 ways
at the same vertex mean and recording the order of the four lobes
(`results/rate_permutations.csv`, one row per assignment):

```text
5040 assignments
  431 give the clinical order (8.6 percent)
  of the 2520 with the occipital rate at or above the frontal one, none
  of those with frontal above occipital, 17 percent
  of those with a ratio of three or more, 75 percent
  the assignment of the reference has the largest possible ratio, 3.3046
```

So the necessary feature is that the anterior groups convert faster than the
posterior ones; the rest of the reference's ranking is not used, the
temporal lobe being first because of the seed, as the uniform run shows.

That condition is a property of this model on this graph and not of tau.
Two limitations are recorded with it. A rate of accumulation measured by
PET is not the coefficient alpha of the equation: it is the observed change
of a burden that transport, local conversion and the concentration already
present produce together, which is why Corti et al. estimate alpha through
an inverse problem rather than reading it from the images. And the regional
pattern of tau accumulation is not stationary: longitudinal tau PET reports
the fastest accumulation in the mesial and inferior temporal cortex before
symptoms and in temporoparietal cortex once they appear (Jack et al., Brain
141 (2018) 1517-1528; Krishnadas, Dore, Robertson et al., eBioMedicine 88
(2023) 104450), so an anterior-to-posterior gradient of conversion is not a
documented property of the protein. The experiment is therefore a
demonstration of the mechanism and of the method, not a calibration for
tau.
The study integrates the same lumped-mass system with Heun's method at a
small step instead of backward Euler at 0.4 years, so its crossing times run
about two years earlier than the stored runs; only the order is used. The figure of this variant is
`results/seeding_patterns_regional_expected.{png,pdf}`, built from the same
strips and the same three-stage rule as the uniform one and collected into
`report/images` beside it: the uniform figure shows what the connectivity
does alone and the regional one what the added biology restores.

## What the amyloid row demonstrates

The amyloid seeding fills the first phase of the expected progression at
`t = 0`, so only the ordering after it can be a result. To see how much the
initial condition carries, `scripts/study-amyloid-single-seed.py` seeds each
of the 54 neocortical vertices alone and repeats the run with the uniform
rate and with the regional field, ordering the four phases (neocortex,
allocortex with the insula and the limbic belt, deep nuclei, brainstem) by
the mean first crossing of `c = 0.5` over their vertices
(`results/amyloid_single_seed.csv`, 108 runs):

```text
                mean phase crossing over the 54 seeds [years]
uniform    neocortex 26.1  allocortex 25.4  deep 24.5  brainstem 25.4
           -> deep < brainstem < allocortex < neocortex
           2.13 of the 6 pairwise orderings right on average, 1 seed of 54
           with all six
regional   neocortex 26.6  allocortex 24.4  deep 25.4  brainstem 27.0
           -> allocortex < deep < neocortex < brainstem
           3.20 of 6 on average, none with all six
```

From a single cortical seed the deep nuclei are reached first, because they
are the hubs of the graph, so neither field reproduces the phases: the
descent of the report figure comes from the neocortex being seeded, not from
the dynamics. The comparison also settles what the regional field does to
this protein. Where the target is a simultaneous onset, as in the seeded
mantle, a uniform rate is the only field that can produce it and any
heterogeneity necessarily breaks it; where the target is an ordering, as
here, the regional field is better than the uniform one (3.20 against 2.13
of six, and 35, 40 and 49 of the 54 seeds for the three orderings after the
cortex, against 16, 35 and 29). The apparent damage to the amyloid row is
therefore a property of the comparison, not of the field.

## The amyloid stations, region by region

The staged renders show the descent only as the centre of the glass brain
lagging behind the mantle, because the deep structures sit inside it in the
lateral view. `results/amyloid_phase_crossings` therefore reports it region
by region: one dot per region at its first crossing of `c = 0.5`, one row
per station, the black bar at the group mean the report quotes. The stations
are composed from the seven groups of benchmark 21: the 54 seeded vertices,
the insula, the limbic belt (which contains the allocortex: entorhinal,
parahippocampal, hippocampi), the subcortical nuclei without the brainstem
and the brainstem itself. The 54 seeded dots collapse onto 4.4 years, which
is the block onset of the mantle; limbic and subcortical spreads overlap
(10.4 to 14.0 and 10.8 to 14.8 years) while their means order as the text
states; the brainstem is last at 17.2.

## Expected-against-obtained composite

The computed rows are rendered from the left side of the head, the frontal
pole to the left of the image, which is the orientation of the clinical
strips above them (medial view, cerebellum to the right); an earlier version
rendered them from the right and was mirrored with respect to the strips.
The verification checks the orientation of every sagittal view of the
report.


The report opens with `results/seeding_patterns_expected.pdf`, the tau row
alone (`--only tau`): the report motivates in the text why amyloid-beta is
not studied on this graph (extracellular spreading against the axonal
pathways the connectome encodes, the choice of Fornari et al. themselves),
while `results/seeding_patterns_full_expected.pdf` keeps both rows as the
record behind that motivation. The figure places the computed row below the
corresponding clinical staging strip of figure 1 of Weickenmeier et al.
(p. 266; the drawings are adopted there from Jucker and Walker). The strips in `reference/` are not redrawn: they are cut from
the article PDF by `scripts/extract-weickenmeier-staging.py`, which pulls
out the single image embedded in the page, a 300 ppi JPEG of about 370
pixels per drawing, without re-rendering it (that resolution is the ceiling
of the published figure), detects the four cartoon rows as bands of
saturated pixels, keeps the two Alzheimer rows after checking their dominant hue, and
drops the label text at the last unsaturated gap. The reference rows are the
literature's expectation and carry no result of ours; every quantitative
element of the composite, the renders and the stage times, comes from the
stored `*_profiles.csv` of this benchmark. The plain
`results/seeding_patterns.pdf` shows the computed rows alone. Sphere radius
and blue intensity both grow with the square root of the concentration, a
mapping stated in the caption, so that the seeded cortex at a tenth of its
final value is visible at the first stage of the full variant's amyloid
row; a region at zero stays a point.

## Reproduce

Run the commands in `commands.txt` from the project root.
