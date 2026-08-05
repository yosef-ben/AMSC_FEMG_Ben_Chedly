#!/usr/bin/env python3

"""Anatomical progression of the misfolded-protein concentration.

Companion to the delayed-conversion figure of Fornari et al. and to the
whole-brain progression figures of Weickenmeier et al. Those works render a
continuum field on a three-dimensional brain mesh; here the substrate is the
metric graph itself, so the concentration is shown on the vertices and edges of
the connectome projected onto an anatomical plane. No volumetric interpolation
is performed and none is implied.

The first two rows vary the conversion rate at the connectivity scaling of the
paper and reproduce its delayed-conversion effect. The third row lowers the
scaling, which is the only way a spatial front becomes visible at all: at unit
scaling the connectome homogenises long before the reaction develops, so every
vertex carries the same concentration at every instant.
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import (load_edges, load_nodes, projection,
                              style_anatomical_axis)

COLOURMAP = plt.cm.inferno
ALPHA_ROOT = Path("output/fisher_kolmogorov/fornari83_alpha")
SCALING_ROOT = Path("output/fisher_kolmogorov/diffusion_scaling/runs")

CASES = (
    (ALPHA_ROOT / "alpha_0p5", r"$\alpha=0.5$" "\n" r"$\rho=1$"),
    (ALPHA_ROOT / "alpha_0p3", r"$\alpha=0.3$" "\n" r"$\rho=1$"),
    (SCALING_ROOT / "rho_0.02", r"$\alpha=0.5$" "\n" r"$\rho=0.02$"),
)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--times", type=float, nargs="+",
                        default=[6.0, 10.0, 14.0, 18.0])
    parser.add_argument("--view", default="sagittal")
    return parser.parse_args()


def read_profiles(path, n_nodes):
    times, states = [], []
    with open(path, newline="") as stream:
        for row in csv.DictReader(stream):
            times.append(float(row["time"]))
            states.append([float(row[f"node_{i}"]) for i in range(n_nodes)])
    return np.array(times), np.array(states)


def main():
    args = arguments()
    nodes = load_nodes()
    edges = load_edges()
    coords = np.array([node["coords"] for node in nodes])
    horizontal, vertical, _, _ = projection(coords, args.view)
    maximum_weight = max(weight for _, _, weight in edges)

    plt.rcParams.update({"font.size": 9.5})
    figure, axes = plt.subplots(len(CASES), len(args.times),
                                figsize=(2.95 * len(args.times),
                                         2.15 * len(CASES)))
    axes = np.atleast_2d(axes)

    for row, (directory, label) in enumerate(CASES):
        # The nodal reference is used throughout: it is the model for which
        # this figure exists in the reference works, and unlike the P1 FEM it
        # stays bounded at every connectivity scaling shown here.
        times, states = read_profiles(directory / "nodal_profiles.csv",
                                      len(nodes))
        for column, instant in enumerate(args.times):
            axis = axes[row, column]
            index = int(np.argmin(np.abs(times - instant)))
            state = states[index]
            # An edge takes the mean of its endpoints, which is what the P1
            # field with one element per connection actually is.
            for source, target, weight in edges:
                axis.plot([horizontal[source], horizontal[target]],
                          [vertical[source], vertical[target]],
                          color=COLOURMAP(0.5 * (state[source]
                                                 + state[target])),
                          linewidth=0.25 + 1.4 * weight / maximum_weight,
                          alpha=0.75, solid_capstyle="round", zorder=1)
            axis.scatter(horizontal, vertical, c=state, cmap=COLOURMAP,
                         vmin=0.0, vmax=1.0, s=26, linewidths=0.4,
                         edgecolors="0.35", zorder=3)
            style_anatomical_axis(axis)
            if row == 0:
                axis.set_title(f"$t = {times[index]:g}$ years", fontsize=10)
            if column == 0:
                axis.text(-0.02, 0.5, label, transform=axis.transAxes,
                          va="center", ha="right", fontsize=10)

    mappable = plt.cm.ScalarMappable(
        cmap=COLOURMAP, norm=matplotlib.colors.Normalize(0.0, 1.0))
    bar = figure.colorbar(mappable, ax=axes, orientation="horizontal",
                          fraction=0.05, pad=0.035, aspect=50)
    bar.set_label("misfolded protein concentration $c$")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220, bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
