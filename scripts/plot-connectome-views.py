#!/usr/bin/env python3

"""Four anatomical views of the 83-region graph, after Fornari et al.

Companion to the brain-network figure of Fornari et al. (their figure 5):
the same two-by-two arrangement of a sagittal, a coronal, a longitudinal and
an unlabelled oblique view, the first three inside a translucent brain
surface and the oblique showing the graph alone, as the printed panel does;
uniform silver vertices, and straight connections whose colour and thickness
both grow with the fibre number of the connection, over a rainbow bar
labelled only at its ends, centred under the views.

Two departures from the published figure are deliberate and are stated in the
caption and in the benchmark record. The vertex coordinates and the pial
surface are the ones distributed with the public Budapest Reference
Connectome viewer, not the MRI-based brain model of the reference, so the
node layout matches anatomy but not their rendering point for point. And the
reference does not state how fibre number maps to colour and width, so a
choice had to be made here: both follow the square root of the fibre number
between the smallest and the largest value. The distribution is strongly
skewed, so the strongest short association bundles keep the warm colours
while the square root keeps the mid-range distinguishable; the printed
figure shows the same broad spread of warm edges, which a purely linear
ramp of so skewed a distribution could not produce.

The fibre numbers are read straight from the edge table of the reconstructed
graph; nothing is rescaled, clipped or reordered.
"""

import argparse
import csv
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import DATA, load_nodes
from render_connectome import (common_box, common_scale, lookup_table,
                               render, show_render)

# Panel order of the published figure: the label position mirrors it, with
# the oblique view unlabelled and the colourbar in its lower-right corner.
# view, label, label side, surface opacity: the oblique panel shows the
# graph alone, as in the printed figure, so its surface is fully transparent
# (the invisible actor still frames the camera like the other panels).
PANELS = (
    ("sagittal_right", "sagittal", "left", 0.08),
    ("coronal", "coronal", "right", 0.08),
    ("longitudinal", "longitudinal", "left", 0.08),
    ("oblique", None, None, 0.0),
)
NODE_SILVER = 0.80
EDGE_RADIUS = (0.16, 1.95)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_fibre_edges():
    """Return (source, target, fibre_number) for the 1130 connections."""
    with open(DATA / "edges.csv", newline="") as stream:
        return [(int(row["source"]), int(row["target"]),
                 float(row["fibre_number"]))
                for row in csv.DictReader(stream)]


def main():
    args = arguments()
    nodes = load_nodes()
    edges = load_fibre_edges()
    coords = np.array([node["coords"] for node in nodes])
    fibres = np.array([fibre for _, _, fibre in edges])
    norm = matplotlib.colors.PowerNorm(0.5, fibres.min(), fibres.max())
    level = [float(norm(fibre)) for fibre in fibres]
    radii = [EDGE_RADIUS[0] + (EDGE_RADIUS[1] - EDGE_RADIUS[0]) * value
             for value in level]

    rainbow = plt.get_cmap("jet")
    edge_table = lookup_table(rainbow, 0.0, 1.0)
    grey_table = lookup_table(plt.cm.gray, 0.0, 1.0)
    scale = common_scale()

    plt.rcParams.update({"font.family": "serif",
                         "font.serif": ["DejaVu Serif"]})
    with tempfile.TemporaryDirectory() as scratch:
        panels = {}
        for view, _, _, opacity in PANELS:
            panels[view] = render(
                Path(scratch) / f"{view}.png", view, coords,
                [NODE_SILVER] * len(nodes), grey_table, node_radius=4.0,
                edges=edges, edge_values=level, edge_table=edge_table,
                edge_radii=radii, scale=scale, surface_opacity=opacity)[0]

        box = common_box(list(panels.values()))
        figure, axes = plt.subplots(2, 2, figsize=(8.6, 7.0))
        for axis, (view, label, side, _) in zip(axes.flat, PANELS):
            show_render(axis, panels[view], box)
            if label:
                axis.text(0.02 if side == "left" else 0.98, 0.015, label,
                          transform=axis.transAxes, fontsize=11.5,
                          ha=side, va="bottom", color="black")

    # The bar of the reference: title above, min and max beside the ends,
    # no numbers. It sits centred under the four views.
    bar_axis = figure.add_axes([0.395, 0.040, 0.21, 0.018])
    matplotlib.colorbar.ColorbarBase(bar_axis, cmap=rainbow,
                                     orientation="horizontal")
    bar_axis.set_xticks([])
    bar_axis.set_title(r"fibre number $n_{IJ}$", fontsize=11.5, pad=5)
    bar_axis.text(-0.06, 0.5, "min", transform=bar_axis.transAxes,
                  fontsize=10.5, ha="right", va="center")
    bar_axis.text(1.06, 0.5, "max", transform=bar_axis.transAxes,
                  fontsize=10.5, ha="left", va="center")

    figure.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.095,
                           hspace=0.06, wspace=0.02)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_dir / "connectome_views.pdf",
                   facecolor="white")
    figure.savefig(args.output_dir / "connectome_views.png", dpi=260,
                   facecolor="white")
    plt.close(figure)

    print(f"Saved connectome_views to {args.output_dir}")
    print(f"  vertices: {len(nodes)}  connections: {len(edges)}")
    print(f"  fibre number range: [{fibres.min():g}, {fibres.max():g}], "
          f"colour and width follow the square-root ramp")
    print(f"  edge radii: [{min(radii):.3f}, {max(radii):.3f}]")


if __name__ == "__main__":
    main()
