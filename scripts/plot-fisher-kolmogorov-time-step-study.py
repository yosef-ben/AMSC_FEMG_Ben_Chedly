#!/usr/bin/env python3

"""Time-step sensitivity of the one-dimensional Fisher-Kolmogorov test.

A sensitivity study, not a convergence study. Weickenmeier et al. report, in
words and without showing it, that the step was varied between dt/4 and 4 dt
around dt = 0.1 and that larger steps produce a spurious increase in spreading.
We repeat the experiment and measure it with the position of the c = 0.5 front.
The finest run is a numerical reference, never an exact solution, so no quantity
here is an error and no slope is drawn.

Panel (b) shows three dimensionless differences from the finest run, all
defined on the profiles at the common final time T:

    e_inf(dt) = max_x |c_dt - c_ref|
    e_2(dt)   = ( (1/|Omega|) * integral_Omega (c_dt - c_ref)^2 dx )^(1/2)
    e_f(dt)   = |x_f(dt) - x_f(ref)|

with Omega = (-1,1). e_2 is a root-mean-square, that is the stored L2 column
divided by sqrt(|Omega|), so that it is dimensionless and comparable with
e_inf; e_f is a displacement in a domain of half-width one, hence also
dimensionless. e_f exists only where the profile still crosses c = 0.5, which
is dt <= 0.1, and the curve stops there.

Styling follows the line figures of the reference: full frame, no grid, sparse
inward ticks, bold labels, the axis name placed after the last tick, direct
coloured series labels and no legend box. No caption is drawn into the image.

Reads time_step_profiles.csv (dt,x,c) and time_step_study.csv. Nothing is
smoothed, fitted or resampled.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LEVEL = 0.5
DOMAIN = 2.0
STEP_COLOUR = {
    0.025: "#12305F", 0.05: "#1F6FB4", 0.1: "#2E9AA6",
    0.2: "#E08214", 0.3: "#CC3311", 0.4: "#7E1212",
}
MEASURE_COLOUR = {"max": "#8A1A1A", "rms": "#2B3A8F", "front": "#137A3C"}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    return parser.parse_args()


def read_profiles(path):
    grouped = defaultdict(list)
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            grouped[float(row["dt"])].append((float(row["x"]),
                                              float(row["c"])))
    return {step: np.array(sorted(points))
            for step, points in grouped.items()}


def right_front(profile, level=LEVEL):
    """Rightmost crossing of `level`, or None if the profile never crosses it."""
    x, c = profile[:, 0], profile[:, 1]
    for k in range(len(c) - 1, 0, -1):
        if (c[k - 1] - level) * (c[k] - level) < 0:
            return x[k - 1] + (level - c[k - 1]) / (c[k] - c[k - 1]) * (
                x[k] - x[k - 1])
    return None


def main():
    args = arguments()
    profiles = read_profiles(args.profiles)
    steps = sorted(profiles)
    finest = steps[0]
    reference = profiles[finest]
    fronts = {step: right_front(profiles[step]) for step in steps}

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 10, "font.weight": "bold",
        "axes.labelweight": "bold", "axes.linewidth": 1.5,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "text.color": "black", "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black",
    })
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.05))

    # ---- (a) the six final profiles ------------------------------------
    axis = axes[0]
    axis.axhline(LEVEL, color="0.6", linewidth=0.9, linestyle=(0, (4, 3)),
                 zorder=1)
    for step in steps:
        profile = profiles[step]
        dashed = step == 0.4
        axis.plot(profile[:, 0], profile[:, 1], color=STEP_COLOUR[step],
                  linewidth=1.9, linestyle=(0, (2.6, 2.2)) if dashed else "-",
                  solid_capstyle="round", zorder=4 if dashed else 3)
    for step in steps:
        if fronts[step] is not None:
            axis.plot(fronts[step], LEVEL, marker="o", markersize=5.5,
                      color=STEP_COLOUR[step], markeredgecolor="white",
                      markeredgewidth=1.0, zorder=5)

    # Colour to time step, and nothing else: the numbers live in the caption.
    # The block sits below the c = 0.5 line, where no curve and no annotation
    # reaches.
    for offset, step in enumerate(steps):
        axis.text(-0.42, 0.40 - 0.074 * offset, rf"$\Delta t = {step:g}$",
                  fontsize=9, fontweight="bold", color=STEP_COLOUR[step],
                  va="center", zorder=6)
    axis.annotate("", xy=(0.60, LEVEL), xytext=(0.26, LEVEL),
                  arrowprops={"arrowstyle": "-|>", "color": "0.35",
                              "linewidth": 1.2, "mutation_scale": 11},
                  zorder=6)
    axis.text(0.0, LEVEL + 0.075, r"increasing $\Delta t$", fontsize=8.5,
              fontweight="bold", color="0.35", ha="center", va="bottom")

    axis.set_xlim(-1, 1)
    axis.set_ylim(-0.03, 1.06)
    axis.set_xticks([-1, -0.5, 0, 0.5, 1])
    axis.set_xticklabels(["-1", "", "0", "", "1"])
    axis.set_yticks([0, LEVEL, 1])
    axis.set_yticklabels(["0", "0.5", "1"])
    axis.set_ylabel(r"concentration $c$", fontsize=10, labelpad=2)
    # Axis name after the last tick, as in the reference figures.
    axis.text(1.0, -0.155, r"  position $x$", fontsize=10, fontweight="bold",
              ha="left", va="center", transform=axis.get_xaxis_transform(),
              clip_on=False)

    # ---- (b) difference from the finest run ----------------------------
    axis = axes[1]
    filled = [step for step in steps if fronts[step] is None]
    # From the first step without a front to the right edge, exactly.
    right_edge = 0.55
    axis.axvspan(min(filled), right_edge, color="#EFE3D2", zorder=0)

    difference = {
        "max": [np.abs(profiles[s][:, 1] - reference[:, 1]).max()
                for s in steps],
        "rms": [np.sqrt(np.trapezoid((profiles[s][:, 1] - reference[:, 1]) ** 2,
                                     profiles[s][:, 0]) / DOMAIN)
                for s in steps],
        "front": [abs(fronts[s] - fronts[finest])
                  if fronts[s] is not None else None for s in steps],
    }
    for key in ("max", "rms", "front"):
        points = [(s, v) for s, v in zip(steps, difference[key])
                  if v is not None and v > 0]
        axis.plot([s for s, _ in points], [v for _, v in points], "-o",
                  color=MEASURE_COLOUR[key], linewidth=1.9, markersize=5,
                  markeredgecolor="white", markeredgewidth=0.8, zorder=3)
    axis.text(0.047, 1.30, r"$e_\infty$", fontsize=10, fontweight="bold",
              color=MEASURE_COLOUR["max"], va="bottom")
    axis.text(0.047, 0.135, r"$e_2$", fontsize=10, fontweight="bold",
              color=MEASURE_COLOUR["rms"], va="bottom")
    axis.text(0.105, 0.052, r"$e_{\mathrm{f}}$", fontsize=10,
              fontweight="bold", color=MEASURE_COLOUR["front"], va="bottom")
    # e_f simply stops where the profile stops crossing c = 0.5; the shaded
    # band and its label say why, so no censored marker is drawn.
    axis.text(0.33, 0.0135, "no  $c = 0.5$  front", fontsize=8.5,
              fontweight="bold", color="#8A5A16", ha="center", va="bottom")

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(0.021, right_edge)
    axis.set_ylim(0.009, 3.0)
    axis.set_xticks(steps)
    axis.set_xticklabels(["0.025", "0.05", "0.1", "0.2", "", "0.4"])
    axis.tick_params(axis="x", labelsize=8.5)
    axis.minorticks_off()
    axis.set_yticks([0.01, 0.1, 1])
    axis.set_yticklabels(["0.01", "0.1", "1"])
    axis.set_ylabel("difference from the finest run", fontsize=10, labelpad=2)
    axis.text(0.55, -0.155, r"  time step $\Delta t$", fontsize=10,
              fontweight="bold", ha="left", va="center",
              transform=axis.get_xaxis_transform(), clip_on=False)

    for letter, axis in zip("ab", axes):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")

    figure.tight_layout(w_pad=3.2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=400, facecolor="white",
                   bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), facecolor="white",
                   bbox_inches="tight")
    print(f"saved {args.output}  (figure size {figure.get_size_inches()})")
    for step in steps:
        front = fronts[step]
        print(f"  dt = {step:<6g} "
              f"x_f = {front:.4f}" if front is not None else
              f"  dt = {step:<6g} no c = 0.5 crossing, "
              f"min c = {profiles[step][:, 1].min():.4f}")
    print(f"  runs with a front: {sum(f is not None for f in fronts.values())}")
    for key in ("max", "rms", "front"):
        values = ", ".join("--" if v is None else f"{v:.6f}"
                           for v in difference[key])
        print(f"  {key:>5}: {values}")


if __name__ == "__main__":
    main()
