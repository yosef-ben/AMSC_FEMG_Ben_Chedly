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

Both rows solve the nodal network model with `alpha = 0.5`, `dt = 0.4` years
and `rho = 0.005`, the scaling benchmark 23 identifies as reproducing the
lobe separation published by Fornari et al. (`Da` of order one hundred). The
nodal model is used deliberately: at this scaling the consistent-mass
metric-graph FEM leaves the physical range and diverges, which is the
validity boundary documented in benchmark 23, while the nodal solution stays
in `[0,1]` throughout, and the plot script asserts it.

Seeding, at `c0 = 0.1`:

- tau: the two entorhinal cortices, the default of the executable;
- amyloid-beta: the 58 vertices of the four cortical lobes, the `neocortex`
  seeding mode of `test_fisher_kolmogorov_fornari83`.

The three stages are selected by a rule, not by eye: the first instants at
which the network mean reaches 10, 40 and 80 percent.

## Results

```text
tau      stages at 10.4, 14.8, 18.8 years
amyloid  stages at  0.8,  4.4,  8.8 years

amyloid-beta, mean group activation (first crossing of c = 0.5):
temporal 4.46, occipital 4.48, parietal 4.80, frontal 4.86,
insular 5.60, limbic 6.33, subcortical 7.75 years
```

The amyloid ordering is the expected one: the four seeded neocortical lobes
activate first and nearly together, then the insula, the limbic belt, and the
subcortical nuclei last, more than three years behind the cortex. For tau the
first activated vertices are the entorhinal seeds and the spreading reaches
the temporal neighbourhood before the rest of the cortex, as in the Braak
sequence; the ordering of the late lobes differs from the clinical one for
the reason recorded in benchmark 23, the weak connectivity of the frontal
pole.

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
