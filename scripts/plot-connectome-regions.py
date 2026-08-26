#!/usr/bin/env python3

"""Anatomical layout and connectivity of the Budapest-83 graph.

The regions figure has two panels. Panel (a) shows the four cortical lobes of
the staging analysis, one brain per lobe in a level three-quarter view
through a near-white pial surface, in the manner of the lobe brains of
Fornari et al.: the vertices of the lobe in the lobe colour of the line
figures, every other vertex in white so the partition stays visible, and
the connections internal to the lobe in a light tint of the same colour. Panel (b) shows the
resulting metric graph in the oblique view without the pial surface, the
connections coloured by connectivity weight: the abstract object the FEM
computes on, where every connection is one edge of unit metric length.

The connectogram figure keeps the seven anatomical groups of Corti et al.

The anatomical panels are rendered with VTK, the engine ParaView is built on.
"""

import argparse
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import (REGION_COLOUR, REGION_ORDER, load_edges,
                              load_nodes, minimal_colourbar)
from figure_style import LOBE_COLOUR
from lobe_scale import classify
from render_connectome import (common_box, common_scale, lookup_table,
                               render, show_render)

# Fixed 2 x 2 order of panel (a): frontal and parietal on top, temporal and
# occipital below.
LOBE_GRID = (("frontal", "parietal"), ("temporal", "occipital"))
# Corti et al. draw their connectogram with the connections above 5% of the
# strongest one; the same threshold is used here.
CONNECTOGRAM_THRESHOLD = 0.05
# The weight row uses the rainbow ramp of the reference's brain-network
# figure, so the two can be read side by side. The group matrix keeps a
# perceptual ramp instead: it is a heatmap with printed values, not a
# render reproducing a published panel.
WEIGHT_MAP = plt.cm.jet
MATRIX_MAP = plt.cm.viridis


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
    # Corti et al. show the principal connections *between* regions, each
    # drawn as a colour gradient from one group to the other; the same
    # construction is used here. Every arc is the quadratic Bezier that
    # matplotlib's "arc3, rad = 0.28" would draw, sampled into segments whose
    # colours interpolate the two group colours; stronger connections are
    # drawn later, hence on top.
    drawn = [(source, target, weight) for source, target, weight in edges
             if nodes[source]["region"] != nodes[target]["region"]
             and weight >= CONNECTOGRAM_THRESHOLD * maximum]
    samples = np.linspace(0.0, 1.0, 33)[:, None]
    for source, target, weight in sorted(drawn, key=lambda edge: edge[2]):
        start = np.array([np.cos(angle[source]), np.sin(angle[source])])
        end = np.array([np.cos(angle[target]), np.sin(angle[target])])
        chord = end - start
        control = (start + end) / 2 + 0.28 * np.array([-chord[1], chord[0]])
        curve = ((1 - samples) ** 2 * start
                 + 2 * samples * (1 - samples) * control
                 + samples ** 2 * end)
        segments = np.stack([curve[:-1], curve[1:]], axis=1)
        first = np.array(matplotlib.colors.to_rgb(
            REGION_COLOUR[nodes[source]["region"]]))
        second = np.array(matplotlib.colors.to_rgb(
            REGION_COLOUR[nodes[target]["region"]]))
        blend = np.linspace(0.0, 1.0, len(segments))[:, None]
        colours = (1 - blend) * first + blend * second
        axis.add_collection(matplotlib.collections.LineCollection(
            segments, colors=colours,
            linewidths=0.35 + 1.9 * weight / maximum,
            alpha=0.8, capstyle="round", zorder=1))
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
    # Square-root ramp, stated in the caption: the reference declares no
    # mapping, the distribution is strongly skewed, and the square root keeps
    # the mid-range readable while leaving the warm colours to the short
    # association bundles. A linear ramp paints nearly everything at the cold
    # end; the logarithmic reading of the same weights lives in the adjacency
    # panel of the topology figure, where cell-level comparison is the point.
    weight_norm = matplotlib.colors.PowerNorm(0.5, weights.min(), weights.max())
    coords = np.array([node["coords"] for node in nodes])
    scale = common_scale()

    weight_table = lookup_table(WEIGHT_MAP, 0.0, 1.0)
    grey_table = lookup_table(plt.cm.gray, 0.0, 1.0)
    edge_level = [float(weight_norm(weight)) for weight in weights]
    edge_radii = [0.22 + 1.45 * value for value in edge_level]

    # Panel (a): the lobe partition of the staging analysis, from the same
    # classification rule the solver uses.
    lobe_of = [classify(node["name"]) for node in nodes]

    with tempfile.TemporaryDirectory() as scratch:
        panels = {}
        for lobe in LOBE_COLOUR:
            # One three-quarter brain per lobe, in the manner of the lobe
            # brains of Fornari et al.: the vertices of the lobe in the
            # lobe colour of the line figures, every other vertex in white
            # so the partition stays visible, the connections internal to
            # the lobe in a light tint of the same colour with the
            # thickness of the weight ramp, and a near-white surface.
            members = [1.0 if group == lobe else 0.0 for group in lobe_of]
            lobe_table = lookup_table(
                matplotlib.colors.ListedColormap(
                    ["#F4F4F4", LOBE_COLOUR[lobe]]), -0.5, 1.5)
            internal = [position for position, (source, target, _)
                        in enumerate(edges)
                        if lobe_of[source] == lobe == lobe_of[target]]
            base = np.array(matplotlib.colors.to_rgb(LOBE_COLOUR[lobe]))
            tint = tuple(1.0 - 0.35 * (1.0 - base))
            print(f"{lobe}: {int(sum(members))} vertices, "
                  f"{len(internal)} internal connections")
            panels["lobe", lobe] = render(
                Path(scratch) / f"lobe_{lobe}.png", "threequarter", coords,
                members, lobe_table, node_radius=2.7,
                edges=[edges[position] for position in internal],
                edge_values=[1.0] * len(internal),
                edge_table=lookup_table(matplotlib.colors.ListedColormap(
                    [tint]), 0.0, 2.0),
                edge_radii=[0.12 + 0.50 * edge_level[position]
                            for position in internal], scale=scale,
                surface_opacity=0.08,
                surface_colour=(0.96, 0.95, 0.94))[0]
        # Panel (b): the metric graph alone, in the oblique view of the
        # reference's brain-network figure, connections by weight with the
        # uniform silver spheres and no pial surface.
        panels["graph"] = render(
            Path(scratch) / "graph_oblique.png", "oblique", coords,
            [0.84] * len(nodes), grey_table, node_radius=3.8, edges=edges,
            edge_values=edge_level, edge_table=weight_table,
            edge_radii=edge_radii, scale=scale, surface_opacity=0.0)[0]

        lobe_box = common_box([panels["lobe", lobe] for lobe in LOBE_COLOUR])
        figure = plt.figure(figsize=(11.8, 5.6))
        grid = figure.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.9],
                                   hspace=0.0, wspace=0.02)
        first = None
        for row, pair in enumerate(LOBE_GRID):
            for column, lobe in enumerate(pair):
                axis = figure.add_subplot(grid[row, column])
                first = first or axis
                show_render(axis, panels["lobe", lobe], lobe_box)
                axis.set_title(lobe, fontsize=10.5, fontweight="bold",
                               color=LOBE_COLOUR[lobe], pad=2)
        oblique = figure.add_subplot(grid[:, 2])
        show_render(oblique, panels["graph"])

    for letter, axis in (("a", first), ("b", oblique)):
        axis.text(0.0, 1.03, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")
    mappable = plt.cm.ScalarMappable(cmap=WEIGHT_MAP, norm=weight_norm)
    bar = minimal_colourbar(figure, mappable, oblique,
                            r"connectivity weight $A_{IJ}$",
                            low=f"{weights.min():.2g}",
                            high=f"{weights.max():.0f}",
                            fraction=0.05, pad=0.03, aspect=45, shrink=0.55)
    bar.ax.minorticks_off()

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
    image = axes[1].imshow(matrix, cmap=MATRIX_MAP, origin="upper")
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
                         color="0.1" if matrix[i, j] > 0.55 * matrix.max()
                         else "white")
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
