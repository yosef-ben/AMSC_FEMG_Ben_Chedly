#!/usr/bin/env python3

"""The connectome summed over the four biomarker lobes, seen from the seed.

Five groups: the four cortical lobes of Fornari et al., the partition over
which every biomarker curve of the chapter is averaged, and the 25 remaining
regions (insular, limbic, subcortical and the brainstem). One line per pair
of groups, with the total connectivity printed on it and the line width
growing with the square root of the value; the star marks the entorhinal
seed, a temporal region in this partition. All sums come from the stored
edge list; nothing is thresholded or dropped.
"""

import argparse
import csv
from math import sqrt
from pathlib import Path

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style

LOBE_KEYS = {
    "temporal": ("temporal", "bankssts", "entorhinal", "fusiform",
                 "parahippocampal"),
    "frontal": ("frontal", "orbitofrontal", "parsopercularis",
                "parsorbitalis", "parstriangularis", "precentral"),
    "parietal": ("parietal", "postcentral", "precuneus", "supramarginal",
                 "paracentral"),
    "occipital": ("cuneus", "occipital", "lingual", "pericalcarine"),
}
POSITION = {
    "frontal": (0.13, 0.76),
    "parietal": (0.60, 0.86),
    "occipital": (0.90, 0.42),
    "temporal": (0.28, 0.16),
    "other": (0.58, 0.42),
}
# Small manual offsets so no two labels collide.
LABEL_SHIFT = {
    ("frontal", "parietal"): (0.0, 0.03),
    ("frontal", "temporal"): (-0.045, 0.0),
    ("frontal", "occipital"): (-0.10, 0.045),
    ("frontal", "other"): (-0.015, 0.035),
    ("occipital", "parietal"): (0.045, 0.02),
    ("parietal", "temporal"): (-0.11, -0.14),
    ("other", "parietal"): (0.038, 0.0),
    ("occipital", "temporal"): (0.0, -0.03),
    ("occipital", "other"): (0.015, 0.035),
    ("other", "temporal"): (-0.01, -0.035),
}


def classify(name):
    lowered = name.lower()
    for lobe, keys in LOBE_KEYS.items():
        if any(key in lowered for key in keys):
            return lobe
    return "other"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path,
                        default=Path("data/connectome/fornari83/nodes.csv"))
    parser.add_argument("--edges", type=Path,
                        default=Path("data/connectome/fornari83/edges.csv"))
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/24_connectome_topology/results/lobe_connectivity.png"))
    return parser.parse_args()


def main():
    args = arguments()
    with open(args.nodes, newline="") as stream:
        groups = {int(row["node_id"]): classify(row["name"])
                  for row in csv.DictReader(stream)}
    counts = {}
    for group in groups.values():
        counts[group] = counts.get(group, 0) + 1
    coupling = {}
    with open(args.edges, newline="") as stream:
        for edge in csv.DictReader(stream):
            one = groups[int(edge["source"])]
            two = groups[int(edge["target"])]
            if one == two:
                continue
            key = tuple(sorted((one, two)))
            coupling[key] = coupling.get(key, 0.0) \
                + float(edge["connectivity_weight"])

    figure_style.apply()
    figure, axis = plt.subplots(figsize=(6.8, 4.9))
    axis.set_aspect("equal")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    heaviest = max(coupling.values())
    for (one, two), value in sorted(coupling.items(), key=lambda kv: kv[1]):
        x = [POSITION[one][0], POSITION[two][0]]
        y = [POSITION[one][1], POSITION[two][1]]
        axis.plot(x, y, color="0.62", solid_capstyle="round",
                  linewidth=0.7 + 3.3 * sqrt(value / heaviest), zorder=1)
        shift = LABEL_SHIFT[(one, two)]
        axis.text((x[0] + x[1]) / 2 + shift[0], (y[0] + y[1]) / 2 + shift[1],
                  f"{value:.1f}", fontsize=8, fontweight="bold",
                  color="0.25", ha="center", va="center", zorder=4)
    for group, (x, y) in POSITION.items():
        colour = figure_style.LOBE_COLOUR.get(group, "0.45")
        axis.add_patch(plt.Circle((x, y), 0.052, facecolor=colour,
                                  edgecolor="white", linewidth=1.4,
                                  zorder=3))
        below = group in ("temporal", "other", "occipital")
        axis.text(x, y - 0.085 if below else y + 0.085,
                  f"{group} ({counts[group]})", fontsize=9.5,
                  fontweight="bold", color=colour, ha="center",
                  va="top" if below else "bottom", zorder=4)
    seed_x, seed_y = POSITION["temporal"]
    axis.plot([seed_x - 0.012], [seed_y + 0.012], marker="*", color="black",
              markersize=13, zorder=5, linestyle="none")
    axis.text(seed_x - 0.075, seed_y + 0.055, "entorhinal seed", fontsize=8.5,
              fontweight="bold", color="black", ha="right", va="center")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)
    figure.savefig(args.output.with_suffix(".pdf"))
    print(f"Written {args.output} and its PDF")
    for (one, two), value in sorted(coupling.items(),
                                    key=lambda kv: -kv[1]):
        print(f"  {one:9s} - {two:9s} {value:8.2f}")


if __name__ == "__main__":
    main()
