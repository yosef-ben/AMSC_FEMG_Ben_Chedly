"""Shared loading and styling for the connectome figures.

The anatomical region of every vertex is read from the CSV written by
``test_fisher_kolmogorov_corti83``, so the figures and the solver always agree
on the classification instead of duplicating it in Python.
"""

from pathlib import Path
import csv

import numpy as np

DATA = Path("data/connectome/fornari83")
REGION_ASSIGNMENT = Path("benchmarks/21_fisher_kolmogorov_corti83/results"
                         "/reaction_coefficients.csv")

# Seven anatomical groups of Corti et al., ordered as in their table 1. The
# hues were selected so that every pair separates by at least 13.9 in OKLab
# (x100) under protanopia, deuteranopia and tritanopia simulation, and by at
# least 17.3 in normal vision.
REGION_COLOUR = {
    "frontal": "#332288",
    "temporal": "#2CA02C",
    "parietal": "#A6761D",
    "insular": "#CC79A7",
    "limbic": "#661100",
    "occipital": "#88CCEE",
    "subcortical": "#F0E442",
}
REGION_ORDER = tuple(REGION_COLOUR)

# Lobes of Fornari et al. figure 7, with their published hue assignment.
LOBE_COLOUR = {
    "temporal": "#2CA02C",
    "frontal": "#D62728",
    "parietal": "#FF7F0E",
    "occipital": "#1F77B4",
}


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
        edges = [(int(row["source"]), int(row["target"]),
                  float(row["connectivity_weight"]))
                 for row in csv.DictReader(stream)]
    return edges


def short_name(name):
    """Strip the FreeSurfer prefix, keeping the hemisphere marker."""
    text = name.replace("ctx-", "").replace("-", " ")
    for prefix, tag in (("rh ", "R "), ("lh ", "L ")):
        if text.startswith(prefix):
            return tag + text[len(prefix):]
    return text


def projection(coords, view):
    """Anatomical projection of Nx3 coordinates onto a viewing plane."""
    if view == "sagittal":          # anterior to the right, superior upwards
        return coords[:, 1], coords[:, 2], "y (anterior)", "z (superior)"
    if view == "axial":             # looking down on the brain
        return coords[:, 0], coords[:, 1], "x (right)", "y (anterior)"
    if view == "coronal":
        return coords[:, 0], coords[:, 2], "x (right)", "z (superior)"
    raise ValueError(f"unknown view {view}")


def draw_edges(axis, nodes, edges, view, weight_scale=1.0, colour="0.55",
               threshold=0.0, zorder=1):
    """Draw the connections, line width proportional to connectivity."""
    coords = np.array([node["coords"] for node in nodes])
    horizontal, vertical, _, _ = projection(coords, view)
    maximum = max(weight for _, _, weight in edges)
    for source, target, weight in edges:
        if weight < threshold * maximum:
            continue
        axis.plot([horizontal[source], horizontal[target]],
                  [vertical[source], vertical[target]],
                  color=colour, linewidth=weight_scale * weight / maximum,
                  solid_capstyle="round", zorder=zorder)


def style_anatomical_axis(axis, xlabel=None, ylabel=None):
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    if xlabel:
        axis.set_xlabel(xlabel, fontsize=8, color="0.4")
    if ylabel:
        axis.set_ylabel(ylabel, fontsize=8, color="0.4")
