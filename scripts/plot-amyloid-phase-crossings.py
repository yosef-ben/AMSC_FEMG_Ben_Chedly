#!/usr/bin/env python3

"""Region-by-region activation of the amyloid-beta staging run.

The staged renders of benchmark 27 show the seeded mantle darkening and the
centre of the glass brain lagging behind it, but the sequence the text
quotes, neocortex, insula, limbic belt, subcortical nuclei, brainstem, is
computed from the stored profiles and no figure showed it. This one does:
one dot per region at the first year its concentration reaches 50 percent,
on one row per group, with the mean of the group marked by a black bar,
which is exactly the statistic the report quotes. The rows are the composite
groups of the report: the 54 seeded neocortical vertices, the insula, the
limbic belt (which contains the allocortex: the entorhinal and
parahippocampal cortices and the hippocampi), the subcortical nuclei without
the brainstem, and the brainstem itself. Everything is read from the stored
results of benchmark 27 and nothing is smoothed or jittered.
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

ROWS = ("seeded neocortex", "insula", "limbic belt", "subcortical nuclei",
        "brainstem")
COLOUR = {"seeded neocortex": "#1F77B4", "insula": "#2CA02C",
          "limbic belt": "#17A2B8", "subcortical nuclei": "#D62CA8",
          "brainstem": "#555555"}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, default=Path(
        "benchmarks/27_connectome_seeding_patterns/results"
        "/amyloid_profiles.csv"))
    parser.add_argument("--regions", type=Path, default=Path(
        "benchmarks/21_fisher_kolmogorov_corti83/results"
        "/reaction_coefficients.csv"))
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/27_connectome_seeding_patterns/results"
        "/amyloid_phase_crossings.png"))
    return parser.parse_args()


def read(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    args = arguments()
    group = {int(row["node_id"]): row["region"]
             for row in read(args.regions)}
    names = {int(row["node_id"]): row["name"] for row in read(args.regions)}
    rows = read(args.profiles)
    times = [float(row["time"]) for row in rows]
    seeded = [k for k in range(83) if float(rows[0][f"node_{k}"]) > 0.0]

    def classify(k):
        if k in seeded:
            return "seeded neocortex"
        if "Brain-Stem" in names[k]:
            return "brainstem"
        if group[k] == "insular":
            return "insula"
        if group[k] == "subcortical":
            return "subcortical nuclei"
        if group[k] == "limbic":
            return "limbic belt"
        return None      # the four unseeded lobar vertices are allocortex,
                         # already counted with the limbic belt by grouping

    crossings = {row: [] for row in ROWS}
    for k in range(83):
        target = classify(k)
        if target is None:
            continue
        series = [float(row[f"node_{k}"]) for row in rows]
        first = next((t for t, v in zip(times, series) if v >= 0.5), None)
        crossings[target].append(first)

    figure_style.apply()
    figure, axis = plt.subplots(figsize=(9.8, 3.1))
    for slot, row in enumerate(ROWS):
        y = len(ROWS) - slot
        values = crossings[row]
        axis.plot(values, [y] * len(values), "o", markersize=6,
                  markerfacecolor=COLOUR[row], markeredgecolor=COLOUR[row],
                  alpha=0.55, linestyle="none", zorder=3)
        mean = sum(values) / len(values)
        axis.plot([mean], [y], marker="|", color="black", markersize=17,
                  markeredgewidth=2.6, linestyle="none", zorder=4)
    slots = {row: len(ROWS) - slot for slot, row in enumerate(ROWS)}
    axis.text(sum(crossings["seeded neocortex"])
              / len(crossings["seeded neocortex"]),
              slots["seeded neocortex"] - 0.42, "group mean", fontsize=8.5,
              fontweight="bold", color="0.35", ha="center")
    axis.set_xlim(0, 19)
    axis.set_xticks([0, 5, 10, 15])
    axis.set_ylim(0.4, len(ROWS) + 0.75)
    axis.set_yticks(list(slots.values()))
    axis.set_yticklabels(list(slots))
    for tick, row in zip(axis.get_yticklabels(), slots):
        tick.set_color(COLOUR[row])
    figure_style.xname(axis, "t [yr]", y=-0.065, fontsize=9)
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)
    figure.savefig(args.output.with_suffix(".pdf"))
    print(f"Written {args.output} and its PDF")
    for row in ROWS:
        values = sorted(crossings[row])
        print(f"  {row:20s} {len(values):2d} regions, {values[0]:5.1f} to "
              f"{values[-1]:5.1f} yr, mean {sum(values) / len(values):6.3f}")


if __name__ == "__main__":
    main()
