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

## Expected-against-obtained composite

The report opens with `results/seeding_patterns_expected.pdf`, which places
each computed row below the corresponding clinical staging strip of figure 1
of Weickenmeier et al. (p. 266; the drawings are adopted there from Jucker
and Walker). The strips in `reference/` are not redrawn: they are cut from
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
final value is visible at the first amyloid stage; a region at zero stays a
point.

## Reproduce

Run the commands in `commands.txt` from the project root.
