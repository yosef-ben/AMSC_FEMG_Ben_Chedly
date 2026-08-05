# Budapest Connectome: 83-Region FreeSurfer Graph

This directory contains the deterministic 83-region graph used as the
reference domain for the Fisher-Kolmogorov application.

It is generated from the Budapest Reference Connectome v3.0 with:

```bash
python3 scripts/prepare-fornari-connectome.py
```

The preprocessing retains fine edges occurring in at least five subjects and
aggregates the 1015 fine parcels by FreeSurfer parent region. Intra-region
edges are removed. For each retained region pair, the fibre counts are summed
and the fibre lengths are averaged for comparison with the published summary
statistics. The edge weight is the sum of the Budapest
`electrical_connectivity_median` values of the fine connections. This is the
dataset's direct estimate of fibre-count-to-length connectivity and correctly
combines parallel fine connections between the same two regions.

The resulting topology and principal statistics reproduce those reported by
Fornari et al.:

```text
fine graph vertices:       1015
fine graph edges:         37477
FreeSurfer regions:          83
region-to-region edges:     1130
unweighted degree range:    6-48
mean fibre number:          40.1619  (paper: 40.2)
mean fibre length:          38.4009 mm (paper: 38.40 mm)
```

`nodes.csv` contains the anatomical coordinates and labels. `edges.csv`
contains the aggregated connectivity data. `summary.json` records the
validation quantities and the source parameters.

This is the canonical application domain. The older files under
`data/connectome/budapest_lcc_*` belong to the preliminary fine-graph
experiment in benchmark 17 and are intentionally kept separate.
