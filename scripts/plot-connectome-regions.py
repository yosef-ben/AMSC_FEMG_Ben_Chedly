#!/usr/bin/env python3

"""Anatomical layout and connectivity of the Budapest-83 graph.

Companion to the graph-discretization figure of Corti et al. and to the brain
network figure of Fornari et al. The MRI and tractography panels of the former
are omitted: neither dataset is available here, and the graph is taken from the
public Budapest Reference Connectome rather than reconstructed from images.

The anatomical panels are rendered with VTK, the engine ParaView is built on.
"""

import argparse
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import (REGION_COLOUR, REGION_ORDER, load_edges,
                              load_nodes, minimal_colourbar)
from render_connectome import (common_box, common_scale, lookup_table,
                               render, show_render)

VIEW_ORDER = ("sagittal", "coronal", "axial")
# Corti et al. draw their connectogram with the connections above 5% of the
# strongest one; the same threshold is used here.
CONNECTOGRAM_THRESHOLD = 0.05
WEIGHT_MAP = plt.cm.magma_r


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


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
    maximum = max(weight for _, _, weight in edges)
    for source, target, weight in edges:
        # Corti et al. show the principal connections *between* regions, so
        # connections internal to a group are left out.
        if nodes[source]["region"] == nodes[target]["region"]:
            continue
        if weight < CONNECTOGRAM_THRESHOLD * maximum:
            continue
        start = np.array([np.cos(angle[source]), np.sin(angle[source])])
        end = np.array([np.cos(angle[target]), np.sin(angle[target])])
        axis.add_patch(FancyArrowPatch(
            start, end, connectionstyle="arc3,rad=0.28", arrowstyle="-",
            linewidth=0.35 + 1.9 * weight / maximum, color="0.45",
            alpha=0.55, zorder=1))
    for node in order:
        axis.plot(np.cos(angle[node]), np.sin(angle[node]), marker="o",
                  markersize=5.5, color=REGION_COLOUR[nodes[node]["region"]],
                  markeredgecolor="0.25", markeredgewidth=0.4, zorder=3)
    for name in REGION_ORDER:
        members = [node for node in order if nodes[node]["region"] == name]
        mean = np.mean([angle[node] for node in members])
        axis.text(1.12 * np.cos(mean), 1.12 * np.sin(mean), name,
                  fontsize=8.5, color=REGION_COLOUR[name],
                  ha="left" if np.cos(mean) > 0.1 else
                     ("right" if np.cos(mean) < -0.1 else "center"),
                  va="bottom" if np.sin(mean) > 0.1 else
                     ("top" if np.sin(mean) < -0.1 else "center"))
    axis.set_xlim(-1.5, 1.5)
    axis.set_ylim(-1.3, 1.3)
    axis.set_aspect("equal")
    axis.set_axis_off()


def main():
    args = arguments()
    nodes = load_nodes()
    edges = load_edges()
    plt.rcParams.update({"font.size": 9.5})

    weights = np.array([weight for _, _, weight in edges])
    weight_norm = matplotlib.colors.LogNorm(weights.min(), weights.max())
    coords = np.array([node["coords"] for node in nodes])
    scale = common_scale()

    # The anatomical groups are categorical, so the vertices carry an index
    # into a lookup table built from the seven colours, not a scalar field.
    region_index = [REGION_ORDER.index(node["region"]) for node in nodes]
    region_table = lookup_table(
        matplotlib.colors.ListedColormap(
            [REGION_COLOUR[name] for name in REGION_ORDER]),
        -0.5, len(REGION_ORDER) - 0.5)
    weight_table = lookup_table(WEIGHT_MAP, 0.0, 1.0)
    grey_table = lookup_table(plt.cm.gray, 0.0, 1.0)
    edge_level = [float(weight_norm(weight)) for weight in weights]
    edge_radii = [0.22 + 1.45 * value for value in edge_level]

    with tempfile.TemporaryDirectory() as scratch:
        panels = {}
        for view in VIEW_ORDER:
            # Top row: vertices by anatomical group, as in Corti et al.
            panels["region", view] = render(
                Path(scratch) / f"region_{view}.png", view, coords,
                region_index, region_table, node_radius=3.6, edges=edges,
                edge_values=[0.66] * len(edges), edge_table=grey_table,
                edge_radii=[0.16] * len(edges), scale=scale)[0]
            # Bottom row: connections by weight, as in Fornari et al.
            panels["weight", view] = render(
                Path(scratch) / f"weight_{view}.png", view, coords,
                [0.80] * len(nodes), grey_table, node_radius=2.6, edges=edges,
                edge_values=edge_level, edge_table=weight_table,
                edge_radii=edge_radii, scale=scale)[0]

        box = common_box(list(panels.values()))
        figure, axes = plt.subplots(2, 3, figsize=(11.8, 6.6))
        for column, view in enumerate(VIEW_ORDER):
            show_render(axes[0, column], panels["region", view], box)
            axes[0, column].set_title(view, fontsize=10, color="0.3", pad=3)
            show_render(axes[1, column], panels["weight", view], box)

    handles = [Line2D([], [], marker="o", linestyle="none", markersize=7,
                      markerfacecolor=REGION_COLOUR[name],
                      markeredgecolor="0.25", markeredgewidth=0.4, label=name)
               for name in REGION_ORDER]
    figure.legend(handles=handles, loc="center", ncol=7, frameon=False,
                  fontsize=8.5, bbox_to_anchor=(0.5, 0.505),
                  handletextpad=0.3, columnspacing=1.1)
    mappable = plt.cm.ScalarMappable(cmap=WEIGHT_MAP, norm=weight_norm)
    bar = minimal_colourbar(figure, mappable, axes[1, :].tolist(),
                            r"connectivity weight $A_{IJ}$",
                            low=f"{weights.min():.2g}",
                            high=f"{weights.max():.0f}",
                            fraction=0.05, pad=0.03, aspect=45, shrink=0.5)
    bar.ax.minorticks_off()
    figure.subplots_adjust(hspace=0.16, wspace=0.02)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_dir / "connectome_regions.pdf",
                   bbox_inches="tight", facecolor="white")
    figure.savefig(args.output_dir / "connectome_regions.png", dpi=260,
                   bbox_inches="tight", facecolor="white")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11.4, 5.0),
                                gridspec_kw={"width_ratios": [1.15, 1.0]})
    connectogram(axes[0], nodes, edges)
    axes[0].set_title("connections above "
                      f"{CONNECTOGRAM_THRESHOLD:.0%} of the strongest",
                      fontsize=10, color="0.3")

    matrix = region_matrix(nodes, edges)
    image = axes[1].imshow(matrix, cmap=WEIGHT_MAP, origin="upper")
    axes[1].set_xticks(range(len(REGION_ORDER)))
    axes[1].set_yticks(range(len(REGION_ORDER)))
    axes[1].set_xticklabels(REGION_ORDER, rotation=45, ha="right")
    axes[1].set_yticklabels(REGION_ORDER)
    axes[1].tick_params(length=0)
    for spine in axes[1].spines.values():
        spine.set_visible(False)
    for i in range(len(REGION_ORDER)):
        for j in range(len(REGION_ORDER)):
            axes[1].text(j, i, f"{matrix[i, j]:.0f}", ha="center",
                         va="center", fontsize=7.5,
                         color="white" if matrix[i, j] > 0.55 * matrix.max()
                         else "0.2")
    axes[1].set_title("connectivity summed over each pair of groups",
                      fontsize=10, color="0.3")
    for tick, name in zip(axes[1].get_xticklabels(), REGION_ORDER):
        tick.set_color(REGION_COLOUR[name])
    for tick, name in zip(axes[1].get_yticklabels(), REGION_ORDER):
        tick.set_color(REGION_COLOUR[name])

    figure.tight_layout()
    figure.savefig(args.output_dir / "connectome_connectogram.pdf",
                   bbox_inches="tight", facecolor="white")
    figure.savefig(args.output_dir / "connectome_connectogram.png", dpi=260,
                   bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"Saved connectome_regions and connectome_connectogram "
          f"to {args.output_dir}")


if __name__ == "__main__":
    main()
