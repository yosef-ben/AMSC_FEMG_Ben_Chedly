#!/usr/bin/env python3

"""Sequential Fisher-Kolmogorov connectome timings, house style.

Two panels: the three measured phases against the unknown count, with the
min-to-max band of the five repetitions, and the per-step cost with its
spread. No asymptotic guide lines are drawn: four points on a graph whose
sparsity pattern changes with the mesh do not identify a complexity law, and
the measured spread is the honest uncertainty to show. Series are labelled
directly in their colour; panel description lives in the LaTeX caption.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style

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

    figure_style.apply()
    figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    axis = axes[0]
    for name, label, colour, marker in PHASES:
        median = column(rows, f"{name}_seconds")
        low = column(rows, f"{name}_min_seconds")
        high = column(rows, f"{name}_max_seconds")
        axis.fill_between(dofs, low, high, color=colour, alpha=0.16,
                          linewidth=0)
        axis.loglog(dofs, median, marker=marker, linestyle="-", color=colour,
                    linewidth=1.8, markersize=5.5,
                    markeredgecolor="white", markeredgewidth=0.7)
        # The total and the time loop nearly coincide, so their end labels
        # are stacked instead of sitting on the curves.
        shift = {"total": 2.6, "solve": 0.42, "assembly": 1.0}[name]
        figure_style.label_series(axis, dofs[-1] * 1.25,
                                  median[-1] * shift, label,
                                  colour, fontsize=9)
    axis.set_xlim(dofs.min() * 0.7, dofs.max() * 9.0)
    axis.set_ylim(2e-4, 8.0)
    axis.set_yticks([1e-3, 1e-2, 1e-1, 1])
    axis.set_yticklabels(["0.001", "0.01", "0.1", "1"])
    axis.minorticks_off()
    axis.set_ylabel("time [s]", labelpad=2)
    figure_style.xname(axis, "degrees of freedom", y=-0.12)

    axis = axes[1]
    low = column(rows, "solve_min_seconds") / 100.0
    high = column(rows, "solve_max_seconds") / 100.0
    axis.errorbar(dofs, per_step,
                  yerr=[per_step - low, high - per_step],
                  fmt="o-", color="#8C2D04", linewidth=1.8, markersize=5.5,
                  markeredgecolor="white", markeredgewidth=0.7,
                  capsize=3, elinewidth=1.1)
    offsets = {"1": (8, -4), "2": (-10, 10), "4": (8, -10), "8": (8, -2)}
    for x, y, row in zip(dofs, per_step, rows):
        cells = row["cells_per_edge"]
        axis.annotate(f"{cells} cell{'' if cells == '1' else 's'}/edge",
                      (x, y), xytext=offsets[cells],
                      textcoords="offset points",
                      fontsize=8, fontweight="bold", color="0.35")
    # The measured inversion between two and four elements per edge: 2.9x
    # the unknowns for a lower cost per step, because the one- and
    # two-element meshes inherit the hub-dominated fill of the connectome.
    plateau = [index for index, row in enumerate(rows)
               if row["cells_per_edge"] in ("2", "4")]
    if len(plateau) == 2:
        first, second = plateau
        axis.annotate("", xy=(dofs[first], per_step[first]),
                      xytext=(dofs[second], per_step[second]),
                      arrowprops={"arrowstyle": "<->", "color": "0.35",
                                  "linewidth": 1.1})
        axis.text(np.sqrt(dofs[first] * dofs[second]) * 1.3,
                  per_step[second] * 0.30,
                  "2.9x the unknowns,\nlower cost per step",
                  fontsize=8.5, fontweight="bold", color="0.35",
                  ha="center", va="top")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(dofs.min() * 0.6, dofs.max() * 3.0)
    axis.set_ylim(per_step.min() * 0.40, per_step.max() * 3.6)
    axis.set_yticks([1e-3, 1e-2])
    axis.set_yticklabels(["0.001", "0.01"])
    axis.minorticks_off()
    axis.set_ylabel("seconds per time step", labelpad=2)
    figure_style.xname(axis, "degrees of freedom", y=-0.12)

    for letter, axis in zip("ab", axes):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")

    figure.tight_layout(w_pad=3.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=240, bbox_inches="tight",
                   facecolor="white")
    if args.output.suffix.lower() != ".pdf":
        figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight",
                       facecolor="white")
    plt.close(figure)
    print(f"Saved {args.output}")
    print(f"  per-step ratio 4-cell / 2-cell: "
          f"{per_step[2] / per_step[1]:.3f}")


if __name__ == "__main__":
    main()
