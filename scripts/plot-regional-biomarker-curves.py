#!/usr/bin/env python3

"""The four lobe biomarkers of the regional-rate tau run, as curves.

Companion panel to the regional staging renders: the same stored run
(entorhinal seed, the seven regional rates rescaled to the vertex mean 0.5,
lumped mass, rho = 0.005, 80 years) read through the biomarker curves of the
four cortical lobes, the presentation of figure 7 of Fornari et al. The
dashed grey curve is the network mean, the horizontal line marks the 50
percent level whose crossings order the lobes and the three vertical lines
mark the instants at which the network mean reaches 10, 40 and 80 percent,
the three stages displayed by the render row above it in the report. Nothing
is computed here: every curve is a column of the stored biomarker table.
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

LOBES = ("temporal", "frontal", "parietal", "occipital")
# Direct series labels, placed on the rising flank of each curve.
LABEL_AT = {
    "temporal": (11.0, 80.0),
    "frontal": (23.2, 62.0),
    "parietal": (33.5, 44.0),
    "occipital": (38.0, 28.0),
}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--biomarkers", type=Path, default=Path(
        "benchmarks/27_connectome_seeding_patterns/results"
        "/tau_regional_biomarkers.csv"))
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/27_connectome_seeding_patterns/results"
        "/regional_biomarker_curves.png"))
    return parser.parse_args()


def main():
    args = arguments()
    with open(args.biomarkers, newline="") as stream:
        rows = list(csv.DictReader(stream))
    times = [float(row["time"]) for row in rows]
    series = {name: [float(row[name]) for row in rows]
              for name in (*LOBES, "global")}

    def first(values, level):
        return next(t for t, v in zip(times, values) if v >= level)

    stages = [first(series["global"], level) for level in (10.0, 40.0, 80.0)]
    crossings = {lobe: first(series[lobe], 50.0) for lobe in LOBES}

    figure_style.apply()
    figure, axis = plt.subplots(figsize=(9.6, 3.1))
    axis.axhline(50.0, color="0.6", linewidth=1.1, linestyle=(0, (4, 3)),
                 zorder=1)
    for stage in stages:
        axis.axvline(stage, color="0.75", linewidth=1.0,
                     linestyle=(0, (2, 3)), zorder=1)
    axis.plot(times, series["global"], color="0.45", linewidth=1.4,
              linestyle="--", zorder=2)
    axis.text(46.0, 88.0, "network mean", fontsize=9, fontweight="bold",
              color="0.45", ha="left", va="center")
    for lobe in LOBES:
        axis.plot(times, series[lobe], color=figure_style.LOBE_COLOUR[lobe],
                  linewidth=2.0, zorder=3)
        x, y = LABEL_AT[lobe]
        axis.text(x, y, lobe, fontsize=9.5, fontweight="bold",
                  color=figure_style.LOBE_COLOUR[lobe], ha="left",
                  va="center")
    # The run extends to 80 years; the axis stops where every curve has
    # saturated, so the panel shares the scale of the other figures.
    axis.set_xlim(0.0, 60.0)
    axis.set_ylim(0.0, 100.0)
    axis.set_yticks([0, 25, 50, 75, 100])
    axis.set_ylabel(r"biomarker abnormality [\%]"
                    if plt.rcParams.get("text.usetex")
                    else "biomarker abnormality [%]", labelpad=2)
    figure_style.xname(axis, "time [years]", y=-0.085, fontsize=10)
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=300, bbox_inches="tight",
                   facecolor="white")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight",
                   facecolor="white")
    print(f"Written {args.output} and its PDF")
    print("  stages (network mean 10/40/80):",
          " ".join(f"{s:.1f}" for s in stages))
    for lobe in LOBES:
        print(f"  {lobe:9s} crosses 50 percent at {crossings[lobe]:5.1f} yr")


if __name__ == "__main__":
    main()
