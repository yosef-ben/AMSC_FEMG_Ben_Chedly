# Budapest Connectome: 83-Region FreeSurfer Graph

This directory contains the deterministic 83-region graph used as the
reference domain for the Fisher-Kolmogorov application.

It is generated from the Budapest Reference Connectome v3.0 with:

```bash
python3 scripts/prepare-fornari-connectome.py
```

The preprocessing retains fine edges occurring in at least five subjects and
aggregates the 1015 fine parcels by FreeSurfer parent region. Intra-region
edges are removed. The five-subject threshold is not stated by Fornari et
al.; it is the only integer that reproduces their fine graph, four keeping
40895 edges and six 34718. The weight of a fine connection is its median
fibre count divided by its median fibre length, and parallel connections
between the same two regions add, as conductances in parallel. Although the
reference describes the fibre quantities as cohort means, it is these median
fields of the public file that reproduce every published graph statistic at
printed precision; the mean fields, and the per-edge
`electrical_connectivity_median` (the median of the per-subject ratios rather
than the ratio of the medians), do not. For each region pair the fibre
counts are also summed and the fibre lengths averaged, for comparison with
the published summary statistics.

The resulting topology and statistics reproduce those reported by
Fornari et al.:

```text
fine graph vertices:       1015
fine graph edges:         37477
FreeSurfer regions:          83
region-to-region edges:     1130
unweighted degree range:    6-48
mean fibre number:          40.1619  (paper: 40.2)
mean fibre length:          38.4009 mm (paper: 38.40 mm)
mean adjacency:             1.5702  (paper: 1.57)
adjacency range:            0.0085 - 35.3221  (paper: 0.01 - 35.32)
weighted degree range:      2.0505 - 127.6435 (paper: 2.1 - 127.6)
```

`nodes.csv` contains the anatomical coordinates and labels. `edges.csv`
contains the aggregated connectivity data. `summary.json` records the
validation quantities and the source parameters.

This is the canonical application domain. The older files under
`data/connectome/budapest_lcc_*` belong to the preliminary fine-graph
experiment in benchmark 17 and are intentionally kept separate.
