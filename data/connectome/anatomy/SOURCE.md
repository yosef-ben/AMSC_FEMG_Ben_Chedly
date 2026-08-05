# Anatomical visualization assets

The files `nodes.json`, `lh.pial`, and
`rh.pial` are downloaded from the public viewer of the
Budapest Reference Connectome:

```text
https://pitgroup.org/apps/connectome/render/data/
```

The official viewer stores one pial coordinate for each of 83 parent regions
and displays subdivisions of the same parent region at the same location. The
preparation script reproduces this mapping and combines the two FreeSurfer pial
surfaces into a legacy binary VTK PolyData file for ParaView.

The downloaded and generated binary assets are not versioned. Recreate them
with:

```bash
python3 scripts/prepare-budapest-anatomy.py --download
```

Use of the dataset should follow the terms and citation guidance published on
the Budapest Reference Connectome website.
