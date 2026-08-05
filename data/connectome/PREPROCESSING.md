# Budapest connectome preprocessing

Source files:

```text
data/budapest_connectome_3.0_209_0_median.graphml
data/budapest_connectome_3.0_209_0_median.csv
```

Regenerate the processed files from the project root with:

```bash
python3 scripts/prepare-budapest-connectome.py
```

The script verifies that GraphML and CSV contain the same undirected edges and
weights. It extracts the largest connected component and removes self-loops;
the original node identifiers remain available in the CSV mapping.

Generated files:

```text
budapest_lcc_fem.txt       topology in the current FEMG graph format
budapest_lcc_nodes.csv     local/original IDs and anatomical metadata
budapest_lcc_edges.csv     remapped edges and connectivity weights
budapest_lcc_summary.txt   human-readable validation summary
budapest_lcc_summary.json  machine-readable validation summary
```

Every metric edge currently has unit length and one finite element. Coordinates
form a deterministic circular layout for technical visualization only; they are
not anatomical coordinates. Connectivity weights remain in the edge table and
will become edge-dependent diffusion coefficients in the nonlinear solver.
