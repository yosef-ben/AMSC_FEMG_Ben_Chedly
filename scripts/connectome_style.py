"""Shared loading and rendering for the connectome figures.

The anatomical region of every vertex is read from the CSV written by
``test_fisher_kolmogorov_corti83``, so the figures and the solver always agree
on the classification instead of duplicating it in Python.

The anatomical panels themselves are produced by ``render_connectome``.
"""

from pathlib import Path
import csv

import numpy as np

DATA = Path("data/connectome/fornari83")
REGION_ASSIGNMENT = Path("benchmarks/21_fisher_kolmogorov_corti83/results"
                         "/reaction_coefficients.csv")

# Seven anatomical groups of Corti et al., ordered as in their table 1. The
# hues are the ones of the legend of their figure 2, sampled from the
# published panels, so the anatomical figures can be read directly against
# the reference. The choice trades colour-vision robustness for
# comparability: under Machado severity-1.0 simulation the smallest pairwise
# OKLab (x100) separation is 3.3 (insular against limbic, tritanopia)
# against 16.5 in normal vision; the renders mitigate this with dark marker
# edges and lightness differences.
REGION_COLOUR = {
    "frontal": "#E40000",
    "temporal": "#ED7D31",
    "parietal": "#FFFF00",
    "insular": "#00E503",
    "limbic": "#00DED2",
    "occipital": "#3432B2",
    "subcortical": "#DE27CE",
}
REGION_ORDER = tuple(REGION_COLOUR)



def load_nodes():
    """Return the 83 vertices with anatomical coordinates and Corti region."""
    with open(DATA / "nodes.csv", newline="") as stream:
        nodes = list(csv.DictReader(stream))
    with open(REGION_ASSIGNMENT, newline="") as stream:
        region = {int(row["node_id"]): row["region"]
                  for row in csv.DictReader(stream)}
    for node in nodes:
        node["node_id"] = int(node["node_id"])
        node["region"] = region[node["node_id"]]
        node["coords"] = np.array([float(node["x"]), float(node["y"]),
                                   float(node["z"])])
    return sorted(nodes, key=lambda node: node["node_id"])


def load_edges():
    """Return the 1130 region-to-region connections and their weights."""
    with open(DATA / "edges.csv", newline="") as stream:
        return [(int(row["source"]), int(row["target"]),
                 float(row["connectivity_weight"]))
                for row in csv.DictReader(stream)]


def short_name(name):
    """Strip the FreeSurfer prefix, keeping the hemisphere marker."""
    text = name.replace("ctx-", "").replace("-", " ")
    for prefix, tag in (("rh ", "R "), ("lh ", "L ")):
        if text.startswith(prefix):
            return tag + text[len(prefix):]
    return text


def minimal_colourbar(figure, mappable, axis, label, low="min", high="max",
                      fontsize=8.5, **kwargs):
    """Horizontal bar labelled only at its ends, as in the reference figures."""
    bar = figure.colorbar(mappable, ax=axis, orientation="horizontal",
                          **kwargs)
    bar.set_ticks([mappable.norm.vmin, mappable.norm.vmax])
    bar.set_ticklabels([low, high])
    bar.ax.tick_params(labelsize=fontsize, length=0)
    bar.set_label(label, fontsize=fontsize)
    bar.outline.set_visible(False)
    return bar
