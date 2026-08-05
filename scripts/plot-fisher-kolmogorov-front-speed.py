#!/usr/bin/env python3

"""Front-speed verification of the one-dimensional Fisher-Kolmogorov solver.

The Fisher-KPP travelling front of c_t - d c_xx = alpha c (1-c) propagates at
the asymptotic speed 2 sqrt(d alpha). Comparing the measured front position
against that value is an independent check of the coupled diffusion-reaction
dynamics, which the space-time sensitivity figure cannot provide.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# One hue per diffusion coefficient; checked for colour-vision deficiency.
DIFFUSION_COLOUR = {1e-4: "#1F77B4", 2e-4: "#FF7F0E", 4e-4: "#2CA02C"}
FIT_WINDOW = (6.0, 18.0)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("speeds", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=1.0)
    return parser.parse_args()


def front_position(xs, cs, level=0.5):
    """Rightmost upward-to-downward crossing of `level`, linearly interpolated."""
    order = np.argsort(xs)
    xs, cs = np.asarray(xs)[order], np.asarray(cs)[order]
    for k in range(len(cs) - 1, 0, -1):
        if cs[k] < level <= cs[k - 1]:
            span = cs[k - 1] - cs[k]
            return xs[k - 1] + (cs[k - 1] - level) / span * (xs[k] - xs[k - 1])
    return float("nan")


def main():
    args = arguments()
    samples = defaultdict(lambda: defaultdict(list))
    with args.profiles.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if float(row["alpha"]) != args.alpha:
                continue
            samples[float(row["diffusion"])][float(row["time"])].append(
                (float(row["x"]), float(row["c"])))
    with args.speeds.open(newline="") as stream:
        speeds = [row for row in csv.DictReader(stream)
                  if float(row["alpha"]) == args.alpha]

    plt.rcParams.update({"font.size": 10})
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.3),
                                constrained_layout=True)

    for diffusion in sorted(samples):
        colour = DIFFUSION_COLOUR[diffusion]
        times = np.array(sorted(samples[diffusion]))
        positions = np.array([
            front_position(*zip(*samples[diffusion][time])) for time in times])
        finite = np.isfinite(positions)
        label = rf"$d={diffusion / 1e-4:g}\cdot10^{{-4}}$"
        axes[0].plot(times[finite], positions[finite], color=colour,
                     linewidth=1.8, label=label)

        window = finite & (times >= FIT_WINDOW[0]) & (times <= FIT_WINDOW[1])
        speed = 2.0 * np.sqrt(diffusion * args.alpha)
        anchor = positions[window][0] - speed * times[window][0]
        axes[0].plot(times[window], speed * times[window] + anchor,
                     color=colour, linewidth=1.1, linestyle="--")

    axes[0].plot([], [], color="0.4", linewidth=1.1, linestyle="--",
                 label=r"slope $2\sqrt{d\,\alpha}$")
    axes[0].set(xlabel=r"time $t$", ylabel=r"right front position, $c=0.5$",
                xlim=(0, 20), ylim=(0, 1.02),
                title=rf"Front propagation at $\alpha={args.alpha:g}$")
    axes[0].grid(True, linewidth=0.4, alpha=0.45)
    axes[0].legend(fontsize=8.5, loc="upper left")

    measured = np.array([float(row["measured_speed"]) for row in speeds])
    theoretical = np.array([float(row["theoretical_speed"]) for row in speeds])
    diffusions = [float(row["diffusion"]) for row in speeds]
    limit = max(measured.max(), theoretical.max()) * 1.15
    axes[1].plot([0, limit], [0, limit], color="0.5", linewidth=1.1,
                 linestyle="--", label=r"$2\sqrt{d\,\alpha}$ (exact)")
    for x, y, diffusion, row in zip(theoretical, measured, diffusions, speeds):
        axes[1].plot(x, y, "o", color=DIFFUSION_COLOUR[diffusion],
                     markersize=9, zorder=3)
        axes[1].annotate(rf"$d={diffusion / 1e-4:g}\cdot10^{{-4}}$" "\n"
                         f"{100 * float(row['relative_error']):.1f}% high",
                         xy=(x, y), xytext=(6, -14),
                         textcoords="offset points", fontsize=8.5,
                         color=DIFFUSION_COLOUR[diffusion])
    axes[1].set(xlabel=r"asymptotic speed $2\sqrt{d\,\alpha}$",
                ylabel="measured front speed", xlim=(0, limit),
                ylim=(0, limit), title="Measured versus asymptotic speed")
    axes[1].grid(True, linewidth=0.4, alpha=0.45)
    axes[1].legend(fontsize=8.5, loc="upper left")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=300)
    figure.savefig(args.output.with_suffix(".png"), dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
