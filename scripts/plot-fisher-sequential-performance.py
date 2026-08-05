#!/usr/bin/env python3

"""Plot the sequential Fisher-Kolmogorov connectome timings."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PHASES = (
    ("total", "total", "#333333", "o"),
    ("solve", "100 time steps", "#D62728", "s"),
    ("assembly", "assembly", "#1F77B4", "^"),
)


def column(rows, name):
    return np.array([float(row[name]) for row in rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    dofs = column(rows, "n_dofs")
    per_step = column(rows, "seconds_per_step")

    plt.rcParams.update({"font.size": 10.5})
    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))

    # No asymptotic guide lines: four points on a graph whose sparsity pattern
    # changes with the mesh do not identify a complexity law, and the measured
    # spread over the five repetitions is the honest uncertainty to show.
    axis = axes[0]
    for name, label, colour, marker in PHASES:
        median = column(rows, f"{name}_seconds")
        low = column(rows, f"{name}_min_seconds")
        high = column(rows, f"{name}_max_seconds")
        axis.fill_between(dofs, low, high, color=colour, alpha=0.18,
                          linewidth=0)
        axis.loglog(dofs, median, marker=marker, linestyle="-", color=colour,
                    linewidth=1.8, markersize=6, label=label)
    axis.set_xlabel("degrees of freedom")
    axis.set_ylabel("time (s)")
    axis.set_title("Sequential Corti-83 baseline\n"
                   "(median of 5 runs, band = min to max)", fontsize=10)
    axis.grid(True, which="both", linewidth=0.45, alpha=0.4)
    axis.legend(frameon=True, fontsize=8.5, loc="center left")

    axis = axes[1]
    low = column(rows, "solve_min_seconds") / 100.0
    high = column(rows, "solve_max_seconds") / 100.0
    axis.errorbar(dofs, per_step,
                  yerr=[per_step - low, high - per_step],
                  fmt="o-", color="#8C2D04", linewidth=1.8, markersize=6,
                  capsize=3, elinewidth=1.0)
    for x, y, row in zip(dofs, per_step, rows):
        cells = int(row["cells_per_edge"])
        axis.annotate(f"{cells} cell{'' if cells == 1 else 's'}/edge\n"
                      f"{row['matrix_nnz']} nnz",
                      (x, y), xytext=(5, -18), textcoords="offset points",
                      fontsize=8)
    plateau = [index for index, row in enumerate(rows)
               if row["cells_per_edge"] in ("2", "4")]
    if len(plateau) == 2:
        first, second = plateau
        axis.annotate("", xy=(dofs[first], per_step[first]),
                      xytext=(dofs[second], per_step[second]),
                      arrowprops={"arrowstyle": "<->", "color": "0.35",
                                  "linewidth": 1.0})
        axis.text(np.sqrt(dofs[first] * dofs[second]),
                  per_step[first] * 1.45,
                  "2.9x the unknowns\nat the same cost",
                  fontsize=8.5, color="0.3", ha="center")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(dofs.min() * 0.6, dofs.max() * 3.0)
    axis.set_ylim(per_step.min() * 0.45, per_step.max() * 3.2)
    axis.set_xlabel("degrees of freedom")
    axis.set_ylabel("seconds per time step")
    axis.set_title("Nonlinear time-step cost is not a function\n"
                   "of the unknown count alone", fontsize=10)
    axis.grid(True, which="both", linewidth=0.45, alpha=0.4)

    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=240, bbox_inches="tight")
    if args.output.suffix.lower() != ".pdf":
        figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
