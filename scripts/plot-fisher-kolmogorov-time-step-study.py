#!/usr/bin/env python3

"""Time-step sensitivity of the one-dimensional Fisher-Kolmogorov problem.

Reproduces the preliminary time-step study of Weickenmeier et al., who varied
the step between dt/4 and 4 dt around dt = 0.1 and reported a spurious increase
in spreading for the larger steps.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REFERENCE_STEP = 0.1


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    return parser.parse_args()


def main():
    args = arguments()
    profiles = defaultdict(list)
    with args.profiles.open(newline="") as stream:
        for row in csv.DictReader(stream):
            profiles[float(row["dt"])].append(
                (float(row["x"]), float(row["c"])))
    with args.summary.open(newline="") as stream:
        summary = sorted(csv.DictReader(stream),
                         key=lambda row: float(row["dt"]))

    steps = sorted(profiles)
    colours = plt.cm.viridis(
        [index / (len(steps) - 1) for index in range(len(steps))])
    plt.rcParams.update({"font.size": 10})
    figure, axes = plt.subplots(1, 3, figsize=(13.4, 4.2),
                                constrained_layout=True)

    for colour, step in zip(colours, steps):
        profile = sorted(profiles[step])
        axes[0].plot([point[0] for point in profile],
                     [point[1] for point in profile],
                     color=colour, linewidth=1.7,
                     label=rf"$\Delta t={step:g}$")
    axes[0].set(xlabel=r"$x$", ylabel=r"$c(x,T)$", xlim=(-1, 1),
                ylim=(-0.02, 1.02),
                title="Final profiles at the common time $T=19.2$")
    axes[0].grid(True, linewidth=0.4, alpha=0.45)
    axes[0].legend(fontsize=8, ncol=2, loc="lower center")

    step_values = [float(row["dt"]) for row in summary]
    front = [float(row["front_position"]) for row in summary]
    mean = [float(row["mean_concentration"]) for row in summary]

    # The two diagnostics are different physical quantities, a position in the
    # domain and a concentration, so each gets its own axis instead of being
    # overlaid on a shared 0-1 scale where their near-coincidence is an
    # accident of units.
    for panel, values, colour, title, ylabel in (
            (axes[1], front, "#1F77B4", "Right front position",
             r"position where $c=0.5$"),
            (axes[2], mean, "#D62728", "Domain mean concentration",
             r"mean of $c(\cdot,T)$ over $[-1,1]$")):
        panel.plot(step_values, values, "o-", color=colour, linewidth=1.7,
                   markersize=6)
        panel.axhline(values[0], color="0.55", linestyle=":", linewidth=1.1)
        panel.text(0.41, values[0] - 0.012,
                   f"finest run: {values[0]:.4f}",
                   fontsize=8, color="0.35", va="top", ha="right")
        panel.axvline(REFERENCE_STEP, color="0.35", linestyle="--",
                      linewidth=1)
        panel.text(REFERENCE_STEP * 1.18, 0.655,
                   rf"reference $\Delta t={REFERENCE_STEP:g}$",
                   fontsize=8, color="0.35")
        panel.set(xlabel=r"$\Delta t$", ylabel=ylabel, xlim=(0, 0.42),
                  ylim=(0.63, 1.05), title=title)
        panel.grid(True, linewidth=0.4, alpha=0.45)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=300)
    figure.savefig(args.output.with_suffix(".png"), dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
