#!/usr/bin/env python3

"""The connectome at the scale of the four biomarker lobes, and what it predicts.

Panel (a) is the connectome compressed over the four cortical lobes of
Fornari et al., the partition over which every biomarker curve of the chapter
is averaged: one line per pair of lobes with the total connectivity printed
on it and the width growing with the square root of the value, the lobes
arranged as in a left sagittal view of the brain, frontal to the left,
parietal above, occipital to the right and temporal below (the star marking
the entorhinal seed), and a short grey stub at every lobe carrying its total
connectivity to the 25 remaining regions (insular, limbic, subcortical and
the brainstem). Only the frontal-occipital and parietal-temporal lines
cross, at the centre, without sharing a vertex. Panel (b) plots the stored lobe
separations of the three models against the lobe-scale Damkohler number
computed from this compressed graph, with the line at one where the
separation sets in. All sums come from the stored edge list and all spreads
from the stored sweeps; nothing is thresholded, fitted or jittered.
"""

import argparse
import csv
from math import sqrt
from pathlib import Path

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style
from lobe_scale import LobeGraph, damkohler_lobe

FIEDLER = 0.772254
# Lobe separation read from figure 7 of Fornari et al. (temporal about 10
# years, occipital about 15.5 years at alpha = 0.5), drawn with a margin of
# 15 percent as in the record figure of the sweep.
REFERENCE_SPREAD_YEARS = 5.5
BENCH = Path("benchmarks/23_fisher_kolmogorov_diffusion_scaling/results")
POSITION = {
    "frontal": (0.02, 0.50),
    "parietal": (0.50, 0.98),
    "occipital": (0.98, 0.50),
    "temporal": (0.50, 0.02),
}
# Lobe-to-lobe labels: fraction along the first-to-second endpoint and shift.
LABEL_POS = {
    ("frontal", "parietal"): (0.50, -0.034, 0.034),
    ("occipital", "parietal"): (0.50, 0.034, 0.034),
    ("frontal", "occipital"): (0.24, 0.0, 0.034),
    ("frontal", "temporal"): (0.50, -0.034, -0.034),
    ("parietal", "temporal"): (0.70, 0.056, 0.0),
    ("occipital", "temporal"): (0.50, 0.034, -0.034),
}
# Grey stubs towards the 25 remaining regions: direction and label offset.
STUB = {
    "frontal": ((-0.11, 0.0), (0.0, -0.052)),
    "parietal": ((0.0, 0.10), (0.075, 0.0)),
    "occipital": ((0.11, 0.0), (0.0, -0.052)),
    "temporal": ((0.0, -0.11), (0.075, 0.0)),
}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/24_connectome_topology/results/lobe_connectivity.png"))
    return parser.parse_args()


def read(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def draw_graph(axis, graph):
    axis.set_aspect("equal")
    axis.set_xlim(-0.24, 1.24)
    axis.set_ylim(-0.24, 1.18)
    axis.axis("off")
    heaviest = max(graph.coupling.values())

    def width(value):
        return 0.7 + 3.3 * sqrt(value / heaviest)

    for (one, two), value in graph.coupling.items():
        if one == "other" or two == "other":
            continue
        (x0, y0), (x1, y1) = POSITION[one], POSITION[two]
        axis.plot([x0, x1], [y0, y1], color="0.62", solid_capstyle="round",
                  linewidth=width(value), zorder=1)
        t, dx, dy = LABEL_POS[(one, two)]
        axis.text(x0 + t * (x1 - x0) + dx, y0 + t * (y1 - y0) + dy,
                  f"{value:.1f}", fontsize=8, fontweight="bold",
                  color="0.25", ha="center", va="center", zorder=4)
    for lobe, (x, y) in POSITION.items():
        value = graph.coupling[tuple(sorted((lobe, "other")))]
        (dx, dy), (lx, ly) = STUB[lobe]
        axis.plot([x, x + dx], [y, y + dy], color="0.62",
                  solid_capstyle="round", linewidth=width(value), zorder=1)
        axis.plot([x + dx], [y + dy], marker="o", markersize=7,
                  color="0.45", linestyle="none", zorder=2)
        axis.text(x + dx + lx, y + dy + ly, f"{value:.1f}", fontsize=8,
                  fontweight="bold", color="0.25", ha="center",
                  va="center", zorder=4)
    for lobe, (x, y) in POSITION.items():
        colour = figure_style.LOBE_COLOUR[lobe]
        axis.add_patch(Circle((x, y), 0.052, facecolor=colour,
                              edgecolor="white", linewidth=1.4, zorder=3))
    axis.text(0.0, 0.60, f"frontal ({graph.counts['frontal']})",
              fontsize=9.5, fontweight="bold",
              color=figure_style.LOBE_COLOUR["frontal"], ha="right",
              va="bottom")
    axis.text(0.435, 1.06, f"parietal ({graph.counts['parietal']})",
              fontsize=9.5, fontweight="bold",
              color=figure_style.LOBE_COLOUR["parietal"], ha="right",
              va="center")
    axis.text(0.96, 0.60, f"occipital ({graph.counts['occipital']})",
              fontsize=9.5, fontweight="bold",
              color=figure_style.LOBE_COLOUR["occipital"], ha="left",
              va="bottom")
    axis.text(0.565, 0.02, f"temporal ({graph.counts['temporal']})",
              fontsize=9.5, fontweight="bold",
              color=figure_style.LOBE_COLOUR["temporal"], ha="left",
              va="center")
    sx, sy = POSITION["temporal"]
    axis.plot([sx], [sy], marker="*", color="black", markersize=15,
              linestyle="none", zorder=6)
    axis.text(0.435, 0.02, "entorhinal seed", fontsize=8.5,
              fontweight="bold", color="black", ha="right", va="center")
    axis.plot([-0.21], [-0.19], marker="o", markersize=7, color="0.45",
              linestyle="none")
    axis.text(-0.17, -0.19, f"the {graph.counts['other']} remaining regions",
              fontsize=8.5, fontweight="bold", color="0.35", ha="left",
              va="center")


def draw_onset(axis, graph, alpha):
    def to_lobe(da_nominal, model):
        rho = alpha / (da_nominal * FIEDLER)
        return damkohler_lobe(rho, alpha, graph.lobe_rate(model))

    series = {
        "nodal": sorted((to_lobe(float(r["damkohler"]), "nodal"),
                         float(r["lobe_spread_years"]))
                        for r in read(BENCH / "diffusion_scaling.csv")),
        "lumped": sorted((to_lobe(float(r["damkohler"]), "lumped"),
                          float(r["lobe_spread_years"]))
                         for r in read(BENCH / "fem_lumped_sweep.csv")
                         if r["scheme"] == "be_lumped"
                         and r["cells_per_edge"] == "1"),
        "consistent": sorted((to_lobe(float(r["damkohler"]), "consistent"),
                              float(r["lobe_spread_years"]))
                             for r in read(BENCH / "fem_consistent_bounded.csv")),
    }
    axis.axvline(1.0, color="0.55", linewidth=1.1, linestyle=(0, (4, 3)),
                 zorder=1)
    axis.axhspan(REFERENCE_SPREAD_YEARS * 0.85, REFERENCE_SPREAD_YEARS * 1.15,
                 color="0.55", alpha=0.16, zorder=0)
    axis.text(0.034, REFERENCE_SPREAD_YEARS * 1.15 + 0.25,
              "separation of the reference", fontsize=8.5, fontweight="bold",
              color="0.4", ha="left", va="bottom")
    style = {"nodal": ("o", "0.25", "-"),
             "lumped": ("s", "#1F77B4", "--"),
             "consistent": ("D", "#D62728", ":")}
    for model, points in series.items():
        marker, colour, line = style[model]
        axis.plot(*zip(*points), marker=marker, markersize=5, color=colour,
                  linewidth=1.4, linestyle=line, zorder=3)
    axis.set_xscale("log")
    axis.set_xlim(0.03, 900)
    axis.set_ylim(0, 13)
    axis.set_yticks([0, 2, 4, 6, 8, 10, 12])
    axis.set_ylabel("lobe separation [years]", labelpad=2)
    figure_style.xname(axis, r"Da$_{\mathrm{lobe}}$", y=-0.075, fontsize=10)
    figure_style.label_series(axis, 16.0, 3.3, "nodal", "0.25", fontsize=9)
    figure_style.label_series(axis, 120.0, 8.2, "FEM, lumped mass",
                              "#1F77B4", fontsize=9)
    figure_style.label_series(axis, 1.1, 9.2, "FEM, consistent mass",
                              "#D62728", fontsize=9)
    axis.text(0.80, 12.2, r"Da$_{\mathrm{lobe}}$ = 1", fontsize=8.5,
              fontweight="bold", color="0.4", ha="right")
    return series


def main():
    args = arguments()
    graph = LobeGraph()
    figure_style.apply()
    figure = plt.figure(figsize=(10.6, 4.8))
    grid = figure.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.22)
    left = figure.add_subplot(grid[0, 0])
    draw_graph(left, graph)
    left.text(0.0, 1.0, "(a)", transform=left.transAxes, fontsize=10.5,
              fontweight="bold", style="italic", va="bottom")
    right = figure.add_subplot(grid[0, 1])
    series = draw_onset(right, graph, args.alpha)
    right.text(0.0, 1.02, "(b)", transform=right.transAxes, fontsize=10.5,
               fontweight="bold", style="italic", va="bottom")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)
    figure.savefig(args.output.with_suffix(".pdf"))
    print(f"Written {args.output} and its PDF")
    for (one, two), value in sorted(graph.coupling.items(),
                                    key=lambda kv: -kv[1]):
        print(f"  {one:9s} - {two:9s} {value:8.2f}")
    low = REFERENCE_SPREAD_YEARS * 0.85
    for model, points in series.items():
        first = next((d for d, s in points if s >= 1.0), float("nan"))
        before = max((d for d, s in points if s < low), default=float("nan"))
        inside = next((d for d, s in points if s >= low), float("nan"))
        print(f"  {model:10s} lobe rate {graph.lobe_rate(model):.3f}, first "
              f"point above one year at Da_lobe = {first:.2f}, reaches the "
              f"reference band between Da_lobe = {before:.1f} and {inside:.1f}")


if __name__ == "__main__":
    main()
