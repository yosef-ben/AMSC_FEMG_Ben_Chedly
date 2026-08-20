#!/usr/bin/env python3

"""Region-by-region activation of the four lobes in the tau staging run.

One dot per region, placed at the first year its concentration reaches 50
percent, on one row per lobe, with the crossing of the lobe biomarker (the
mean of the group) marked by a black bar. The figure shows why the order of
the lobe means differs from the order of the direct couplings: the mean of a
lobe is set by its late regions, the postcentral and paracentral tail in the
parietal lobe and the two poles in the frontal one. Everything is read from
the stored results of benchmark 27 and nothing is smoothed or jittered.
"""

import argparse
import csv
from pathlib import Path

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style

# The four biomarker lobes of Fornari et al., by FreeSurfer name, as in
# test_fisher_kolmogorov_fornari83.cpp.
LOBE_KEYS = {
    "temporal": ("temporal", "bankssts", "entorhinal", "fusiform",
                 "parahippocampal"),
    "frontal": ("frontal", "orbitofrontal", "parsopercularis",
                "parsorbitalis", "parstriangularis", "precentral"),
    "parietal": ("parietal", "postcentral", "precuneus", "supramarginal",
                 "paracentral"),
    "occipital": ("cuneus", "occipital", "lingual", "pericalcarine"),
}
# Rows from top to bottom, the activation order of the means.
ROWS = ("temporal", "occipital", "parietal", "frontal")


def classify(name):
    lowered = name.lower()
    for lobe, keys in LOBE_KEYS.items():
        if any(key in lowered for key in keys):
            return lobe
    return None


def crossing(times, values, level=50.0):
    for k in range(1, len(values)):
        if values[k - 1] < level <= values[k]:
            span = values[k] - values[k - 1]
            return times[k - 1] + (level - values[k - 1]) / span * (
                times[k] - times[k - 1])
    return float("nan")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(
        "benchmarks/27_connectome_seeding_patterns/results"))
    parser.add_argument("--nodes", type=Path,
                        default=Path("data/connectome/fornari83/nodes.csv"))
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/27_connectome_seeding_patterns/results"
        "/lobe_crossings.png"))
    return parser.parse_args()


def read(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    args = arguments()
    nodes = read(args.nodes)
    profiles = read(args.results / "tau_profiles.csv")
    biomarkers = read(args.results / "tau_biomarkers.csv")
    times = [float(row["time"]) for row in profiles]

    per_lobe = {lobe: [] for lobe in ROWS}
    for node in nodes:
        lobe = classify(node["name"])
        if lobe is None:
            continue
        values = [100.0 * float(row[f"node_{node['node_id']}"])
                  for row in profiles]
        per_lobe[lobe].append((crossing(times, values), node["name"]))
    mean_times = [float(row["time"]) for row in biomarkers]
    means = {lobe: crossing(mean_times,
                            [float(row[lobe]) for row in biomarkers])
             for lobe in ROWS}

    figure_style.apply()
    figure, axis = plt.subplots(figsize=(9.8, 3.3))
    for slot, lobe in enumerate(ROWS):
        y = len(ROWS) - slot
        colour = figure_style.LOBE_COLOUR[lobe]
        xs = [value for value, _ in per_lobe[lobe]]
        axis.plot(xs, [y] * len(xs), "o", markersize=6,
                  markerfacecolor=colour, markeredgecolor=colour,
                  alpha=0.55, linestyle="none", zorder=3)
        axis.plot([means[lobe]], [y], marker="|", color="black",
                  markersize=17, markeredgewidth=2.6, linestyle="none",
                  zorder=4)
    rows = {lobe: len(ROWS) - slot for slot, lobe in enumerate(ROWS)}
    label = dict(fontsize=8.5, fontweight="bold", color="0.35", ha="center")
    axis.text(4.3, rows["temporal"] + 0.27, "entorhinal seed", **label)
    axis.text(means["temporal"], rows["temporal"] - 0.38, "lobe mean",
              **label)
    axis.text(25.5, rows["parietal"] + 0.27, "postcentral, paracentral",
              **label)
    axis.text(30.4, rows["frontal"] + 0.27, "frontal poles", **label)
    axis.set_xlim(0, 33)
    axis.set_xticks([0, 5, 10, 15, 20, 25, 30])
    axis.set_ylim(0.35, 4.78)
    axis.set_yticks(list(rows.values()))
    axis.set_yticklabels(list(rows))
    for tick, lobe in zip(axis.get_yticklabels(), rows):
        tick.set_color(figure_style.LOBE_COLOUR[lobe])
    figure_style.xname(axis, "t [yr]", y=-0.06, fontsize=9)
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)
    figure.savefig(args.output.with_suffix(".pdf"))
    print(f"Written {args.output} and its PDF")
    for lobe in ROWS:
        xs = sorted(per_lobe[lobe])
        print(f"  {lobe:9s} {len(xs):2d} regions, {xs[0][0]:.2f} to "
              f"{xs[-1][0]:.2f} yr, mean curve crossing {means[lobe]:.2f}")


if __name__ == "__main__":
    main()
