#!/usr/bin/env python3

"""Regional concentrations of the deterministic Corti model.

One curve per anatomical group, in the shared group palette, plus the dashed
network mean; every curve is labelled directly at its right end in its own
colour, following the house conventions of figure_style. Nothing is smoothed
or resampled: the points are the stored samples.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style
from connectome_style import REGION_COLOUR


def arguments():
    parser = argparse.ArgumentParser(
        description="Plot regional Fisher-Kolmogorov concentrations."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("regional_averages.png"))
    return parser.parse_args()


def main():
    args = arguments()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    time = [float(row["time"]) for row in rows]
    names = [
        "frontal", "temporal", "parietal", "insular", "limbic",
        "occipital", "subcortical",
    ]

    figure_style.apply()
    fig, axis = plt.subplots(figsize=(7.2, 4.2))
    finals = {}
    for name in names:
        values = [float(row[name]) for row in rows]
        axis.plot(time, values, color=REGION_COLOUR[name], linewidth=1.9,
                  solid_capstyle="round")
        finals[name] = values[-1]
    mean = [float(row["global"]) for row in rows]
    axis.plot(time, mean, color="black", linewidth=2.0,
              linestyle=(0, (5, 3)))
    finals["network mean"] = mean[-1]

    # Direct labels at the right ends. The insular and subcortical curves end
    # 0.0013 apart, so their labels are nudged symmetrically; every other
    # label sits at its own final value.
    nudge = {"insular": +0.006, "subcortical": -0.006,
             "temporal": +0.002, "network mean": -0.002,
             "limbic": +0.004, "frontal": -0.004}
    for name, value in finals.items():
        colour = "black" if name == "network mean" else REGION_COLOUR[name]
        figure_style.label_series(axis, time[-1] + 0.25,
                                  value + nudge.get(name, 0.0),
                                  name, colour)

    axis.set_xlim(time[0], time[-1])
    axis.set_ylim(0.0, 0.225)
    axis.set_xticks([0, 5, 10, 15, 20])
    axis.set_yticks([0.0, 0.05, 0.10, 0.15, 0.20])
    axis.set_yticklabels(["0", "0.05", "0.10", "0.15", "0.20"])
    axis.set_ylabel("mean concentration", labelpad=2)
    figure_style.xname(axis, "time [years]")

    fig.tight_layout()
    fig.subplots_adjust(right=0.82)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=260, facecolor="white",
                bbox_inches="tight")
    if args.output.suffix != ".pdf":
        fig.savefig(args.output.with_suffix(".pdf"), facecolor="white",
                    bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output}")
    for name in names + ["network mean"]:
        print(f"  {name:14s} final {finals[name]:.6f}")


if __name__ == "__main__":
    main()
