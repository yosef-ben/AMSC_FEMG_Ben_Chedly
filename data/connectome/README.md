# Budapest connectome preprocessing

This folder contains the deterministic preprocessing output used to turn the
Budapest Reference Connectome into an input domain for FEMG.

Run:

```bash
python3 scripts/prepare-budapest-connectome.py
python3 scripts/prepare-budapest-anatomy.py --download
```

The first command validates the GraphML and CSV files, extracts their largest
connected component, and writes the FEM topology and metadata. The second maps
each retained fine node to the pial coordinate of its Budapest parent region
and converts the two FreeSurfer pial surfaces into
`anatomy/brain_surface.vtk`.

Metric edge lengths belong to the mathematical model. Anatomical coordinates
are used only by the Fisher-Kolmogorov VTK exporter and never alter the FEM
matrices. The downloaded viewer assets are intentionally ignored by Git and can
be recreated with the second command.
