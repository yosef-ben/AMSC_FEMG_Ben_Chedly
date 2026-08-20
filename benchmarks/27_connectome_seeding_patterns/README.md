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
- amyloid-beta: the 58 vertices of the four cortical lobes, the `neocortex`
  seeding mode of `test_fisher_kolmogorov_fornari83`.

The three stages are selected by a rule, not by eye: the first instants at
which the network mean reaches 10, 40 and 80 percent.

## Results

```text
tau      stages at 15.2, 21.6, 27.6 years
amyloid  stages at  0.8,  4.8, 11.2 years

amyloid-beta, mean group activation (first crossing of c = 0.5):
temporal 4.4, parietal 4.4, occipital 4.4, frontal 4.4,
insular 9.6, limbic 9.78, subcortical 12.98 years

tau, lobe 50-percent crossings: temporal 16.50, occipital 21.65,
parietal 22.95, frontal 25.87 years
```

The amyloid ordering is the expected one: the four seeded neocortical lobes
activate first and together, then the insula, the limbic belt, and the
subcortical nuclei last, more than eight years behind the cortex. For tau
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
the article PDF by `scripts/extract-weickenmeier-staging.py`, which renders
the page at 600 dpi, detects the four cartoon rows as bands of saturated
pixels, keeps the two Alzheimer rows after checking their dominant hue, and
drops the label text at the last unsaturated gap. The reference rows are the
literature's expectation and carry no result of ours; every quantitative
element of the composite, the renders and the stage times, comes from the
stored `*_profiles.csv` of this benchmark. The plain
`results/seeding_patterns.pdf` shows the computed rows alone.

## Reproduce

Run the commands in `commands.txt` from the project root.
