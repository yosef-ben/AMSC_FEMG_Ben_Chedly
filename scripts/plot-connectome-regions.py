#!/usr/bin/env python3

"""Anatomical layout and regional connectivity of the Budapest-83 graph.

Companion to the graph-discretization figure of Corti et al. Their panels for
the MRI surface and the DWI tractography are omitted: neither dataset is
available here, and the graph is taken from the public Budapest Reference
Connectome rather than reconstructed from images.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import (REGION_COLOUR, REGION_ORDER, draw_edges,
                              load_edges, load_nodes, projection,
                              style_anatomical_axis)

# Corti et al. draw their connectogram with the connections above 5% of the
# strongest one; the same threshold is used here.
CONNECTOGRAM_THRESHOLD = 0.05


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def scatter_nodes(axis, nodes, view, size=42, zorder=3):
    coords = np.array([node["coords"] for node in nodes])
    horizontal, vertical, xlabel, ylabel = projection(coords, view)
    colours = [REGION_COLOUR[node["region"]] for node in nodes]
    axis.scatter(horizontal, vertical, s=size, c=colours, linewidths=0.5,
                 edgecolors="0.25", zorder=zorder)
    style_anatomical_axis(axis, xlabel, ylabel)


def region_matrix(nodes, edges):
    index = {name: position for position, name in enumerate(REGION_ORDER)}
    matrix = np.zeros((len(REGION_ORDER), len(REGION_ORDER)))
    for source, target, weight in edges:
        i = index[nodes[source]["region"]]
        j = index[nodes[target]["region"]]
        matrix[i, j] += weight
        if i != j:
            matrix[j, i] += weight
    return matrix


def connectogram(axis, nodes, edges):
    """Circular layout grouped by region, strongest connections only."""
    order = sorted(range(len(nodes)),
                   key=lambda i: (REGION_ORDER.index(nodes[i]["region"]),
                                  nodes[i]["name"]))
    angle = {node: 2 * np.pi * position / len(order)
             for position, node in enumerate(order)}
    radius = 1.0
    maximum = max(weight for _, _, weight in edges)
    for source, target, weight in edges:
        # Corti et al. show the principal connections *between* regions, so
        # connections internal to a group are left out.
        if nodes[source]["region"] == nodes[target]["region"]:
            continue
        if weight < CONNECTOGRAM_THRESHOLD * maximum:
            continue
        start = np.array([radius * np.cos(angle[source]),
                          radius * np.sin(angle[source])])
        end = np.array([radius * np.cos(angle[target]),
                        radius * np.sin(angle[target])])
        # Bend towards the centre so that long-range connections stay legible.
        axis.add_patch(FancyArrowPatch(
            start, end, connectionstyle="arc3,rad=0.28", arrowstyle="-",
            linewidth=0.35 + 1.9 * weight / maximum, color="0.45",
            alpha=0.55, zorder=1))
    for node in order:
        position = np.array([radius * np.cos(angle[node]),
                             radius * np.sin(angle[node])])
        axis.plot(*position, marker="o", markersize=5.5,
                  color=REGION_COLOUR[nodes[node]["region"]],
                  markeredgecolor="0.25", markeredgewidth=0.4, zorder=3)
    # Name each group once, at the mean angle of its vertices.
    for name in REGION_ORDER:
        members = [node for node in order if nodes[node]["region"] == name]
        if not members:
            continue
        mean = np.mean([angle[node] for node in members])
        axis.text(1.12 * np.cos(mean), 1.12 * np.sin(mean), name,
                  fontsize=8.5, color=REGION_COLOUR[name],
                  ha="left" if np.cos(mean) > 0.1 else
                     ("right" if np.cos(mean) < -0.1 else "center"),
                  va="bottom" if np.sin(mean) > 0.1 else
                     ("top" if np.sin(mean) < -0.1 else "center"))
    axis.set_xlim(-1.5, 1.5)
    axis.set_ylim(-1.3, 1.3)
    style_anatomical_axis(axis)


def main():
    args = arguments()
    nodes = load_nodes()
    edges = load_edges()
    plt.rcParams.update({"font.size": 9.5})

    # Anatomical layout, two views with and without the connections.
    figure, axes = plt.subplots(1, 3, figsize=(12.4, 4.2))
    scatter_nodes(axes[0], nodes, "sagittal")
    for axis, view in ((axes[1], "sagittal"), (axes[2], "axial")):
        draw_edges(axis, nodes, edges, view, weight_scale=2.4, colour="0.65")
        scatter_nodes(axis, nodes, view, size=26)
    # Equal aspect gives the three panels different heights, so the titles are
    # placed at a common height in figure coordinates instead of per axes.
    for axis, label in zip(axes, ("(a) vertices by anatomical group",
                                  "(b) graph, sagittal view",
                                  "(c) graph, axial view")):
        box = axis.get_position()
        figure.text(box.x0 + box.width / 2, 0.95, label, fontsize=10,
                    ha="center", va="bottom")

    handles = [Line2D([], [], marker="o", linestyle="none", markersize=7,
                      markerfacecolor=REGION_COLOUR[name],
                      markeredgecolor="0.25", markeredgewidth=0.4, label=name)
               for name in REGION_ORDER]
    figure.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
                  fontsize=9, bbox_to_anchor=(0.5, -0.01))
    figure.tight_layout(rect=(0, 0.06, 1, 0.94))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_dir / "connectome_regions.pdf",
                   bbox_inches="tight")
    figure.savefig(args.output_dir / "connectome_regions.png", dpi=220,
                   bbox_inches="tight")
    plt.close(figure)

    # Connectogram and aggregated regional connectivity.
    figure, axes = plt.subplots(1, 2, figsize=(11.4, 5.0),
                               gridspec_kw={"width_ratios": [1.15, 1.0]})
    connectogram(axes[0], nodes, edges)
    axes[0].set_title("(a) connections above "
                      f"{CONNECTOGRAM_THRESHOLD:.0%} of the strongest",
                      fontsize=10)

    matrix = region_matrix(nodes, edges)
    image = axes[1].imshow(matrix, cmap="magma_r", origin="upper")
    axes[1].set_xticks(range(len(REGION_ORDER)))
    axes[1].set_yticks(range(len(REGION_ORDER)))
    axes[1].set_xticklabels(REGION_ORDER, rotation=45, ha="right")
    axes[1].set_yticklabels(REGION_ORDER)
    for i in range(len(REGION_ORDER)):
        for j in range(len(REGION_ORDER)):
            axes[1].text(j, i, f"{matrix[i, j]:.0f}", ha="center",
                         va="center", fontsize=7.5,
                         color="white" if matrix[i, j] > 0.55 * matrix.max()
                         else "0.2")
    axes[1].set_title("(b) total connectivity between groups", fontsize=10)
    bar = figure.colorbar(image, ax=axes[1], fraction=0.046, pad=0.03)
    bar.set_label("summed connectivity weight", fontsize=9)
    for tick, name in zip(axes[1].get_xticklabels(), REGION_ORDER):
        tick.set_color(REGION_COLOUR[name])
    for tick, name in zip(axes[1].get_yticklabels(), REGION_ORDER):
        tick.set_color(REGION_COLOUR[name])

    figure.tight_layout()
    figure.savefig(args.output_dir / "connectome_connectogram.pdf",
                   bbox_inches="tight")
    figure.savefig(args.output_dir / "connectome_connectogram.png", dpi=220,
                   bbox_inches="tight")
    plt.close(figure)
    print(f"Saved connectome_regions and connectome_connectogram "
          f"to {args.output_dir}")


if __name__ == "__main__":
    main()
