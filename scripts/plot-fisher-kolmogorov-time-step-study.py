#!/usr/bin/env python3
"""Time-step sensitivity of the one-dimensional Fisher-Kolmogorov test.

Two panels: (a) the six final profiles at T = 19.2 with the c = 0.5 fronts
marked, (b) one convergence panel in which the front position and the two
solution norms are all expressed as a *relative error against the dt = 0.025
run*, so that they legitimately share a single dimensionless axis.

Every mark is read from
  benchmarks/18_fisher_kolmogorov_1d_sensitivity/results/time_step_profiles.csv
  benchmarks/18_fisher_kolmogorov_1d_sensitivity/results/time_step_study.csv
Nothing is smoothed, fitted to the curves, rescaled or resampled.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

RESULTS = None

REFERENCE = 0.025          # finest run, used as the numerical reference
FILLED_FROM = 0.2          # smallest tested step whose solution fills [-1,1]
FIT_PAIR = (0.05, 0.1)     # the only two steps with a genuine front and error

# Steps are coloured dark-cool to bright-warm with increasing dt.  The three
# cool hues are the runs that still carry a front, the three warm ones the runs
# that have filled the domain; the same cut is the shaded band in (b).
COLOUR = {
    0.025: "#0E1E46",
    0.05:  "#1F5FA8",
    0.1:   "#2196A8",
    0.2:   "#E08214",
    0.3:   "#D3402A",
    0.4:   "#7E1212",
}
INK = "#1A1A1A"
BAND = "#F2E4D2"


def load():
    profiles = defaultdict(list)
    with (RESULTS / "time_step_profiles.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            profiles[float(row["dt"])].append((float(row["x"]),
                                               float(row["c"])))
    grids = {}
    for step, points in profiles.items():
        points.sort()
        grids[step] = (np.array([p[0] for p in points]),
                       np.array([p[1] for p in points]))
    with (RESULTS / "time_step_study.csv").open(newline="") as stream:
        summary = {float(row["dt"]): {k: float(v) for k, v in row.items()}
                   for row in csv.DictReader(stream)}
    return grids, summary


def crossing(x, c, level=0.5):
    """Right-hand c = level crossing of the P1 interpolant, or None."""
    for i in range(len(x) - 1, 0, -1):
        a, b = c[i - 1], c[i]
        if (a - level) * (b - level) < 0.0:
            return x[i - 1] + (level - a) / (b - a) * (x[i] - x[i - 1])
    return None


def diagnostics(grids, steps):
    xr, cref = grids[REFERENCE]
    ref_norm = np.sqrt(np.trapezoid(cref ** 2, xr))
    front_ref = crossing(xr, cref)
    rel = {"l2": {}, "inf": {}, "front": {}}
    fronts, minima = {}, {}
    for step in steps:
        x, c = grids[step]
        assert np.allclose(x, xr)
        rel["l2"][step] = np.sqrt(np.trapezoid((c - cref) ** 2, x)) / ref_norm
        rel["inf"][step] = np.max(np.abs(c - cref)) / np.max(np.abs(cref))
        fronts[step] = crossing(x, c)
        minima[step] = c.min()
        # No crossing -> the front has left [-1,1]; the value is then censored
        # at the domain end x = 1, exactly what the solver wrote in the CSV.
        rel["front"][step] = abs((fronts[step] if fronts[step] is not None
                                  else 1.0) - front_ref) / front_ref
    # value of each measure for the trivial filled state c == 1
    one = np.ones_like(cref)
    ceiling = {
        "l2": np.sqrt(np.trapezoid((one - cref) ** 2, xr)) / ref_norm,
        "inf": np.max(np.abs(one - cref)) / np.max(np.abs(cref)),
        "front": abs(1.0 - front_ref) / front_ref,
    }
    order = {k: np.log(rel[k][FIT_PAIR[1]] / rel[k][FIT_PAIR[0]])
             / np.log(FIT_PAIR[1] / FIT_PAIR[0]) for k in rel}
    return rel, fronts, minima, ceiling, order, front_ref


def profiles_panel(axis, grids, steps, fronts, minima):
    axis.axhline(0.5, color="#C9C9C9", linewidth=1.4, linestyle=(0, (5, 4)),
                 zorder=1)
    for step in steps:
        x, c = grids[step]
        # dt = 0.3 and dt = 0.4 differ by less than 4e-5 in c and would be one
        # line; the coarser of the two is dashed so that both stay visible.
        style = dict(linestyle=(0, (2.2, 2.4)), linewidth=2.4) if step == 0.4 \
            else dict(linewidth=3.0, solid_capstyle="round")
        axis.plot(x, c, color=COLOUR[step], zorder=3 + (step == 0.4), **style)

    axis.set_xlim(-1.06, 1.06)
    axis.set_ylim(-0.045, 1.06)
    axis.set_xticks([-1, -0.5, 0, 0.5, 1])
    axis.set_yticks([0, 0.5, 1])
    axis.set_xticklabels(["\u22121", "\u22120.5", "0", "0.5", "1"])
    axis.set_yticklabels(["0", "0.5", "1"])
    axis.set_xlabel("$x$", fontsize=15, fontweight="bold", labelpad=1)
    axis.set_ylabel("concentration $c(x, T)$", fontsize=13.5,
                    fontweight="bold", labelpad=2)
    axis.tick_params(labelsize=13)

    # the c = 0.5 fronts, and the direction in which they move
    axis.annotate("", xy=(0.925, 0.5), xytext=(fronts[REFERENCE], 0.5),
                  arrowprops=dict(arrowstyle="-|>", linewidth=1.6,
                                  color="#8C8C8C", shrinkA=0, shrinkB=0),
                  zorder=2)
    axis.text(0.945, 0.5, r"$\Delta t\,\uparrow$", fontsize=13.5,
              color="#6E6E6E", va="center", ha="left", fontweight="bold")
    for step in (0.025, 0.05, 0.1):
        axis.plot([fronts[step]], [0.5], "o", color=COLOUR[step],
                  markersize=10, markeredgecolor="white", markeredgewidth=1.6,
                  zorder=6)
    axis.text(-1.04, 0.535, r"$c = 0.5$", fontsize=10.5, color="#8C8C8C",
              ha="left", va="bottom", style="italic")

    rows = [
        (0.025, r"front $x_{\rm f} = 0.6889$   (reference)"),
        (0.05, r"front $x_{\rm f} = 0.7268$"),
        (0.1, r"front $x_{\rm f} = 0.8145$"),
        (0.2, r"no front,  $c_{\min} = 0.9320$"),
        (0.3, r"no front,  $c_{\min} = 0.99996$"),
        (0.4, r"no front,  $c_{\min} = 0.99999994$"),
    ]
    for index, (step, tail) in enumerate(rows):
        axis.text(-0.63, 0.795 - 0.125 * index,
                  rf"$\Delta t = {step:g}$     {tail}",
                  fontsize=10.5, color=COLOUR[step], fontweight="bold",
                  ha="left", va="center")
    axis.text(-0.63, 0.795 + 0.115,
              r"$x_{\rm f}$ is where the profile crosses $c = 0.5$",
              fontsize=10.5, color="#6E6E6E", ha="left", va="center",
              style="italic")
    # how closely the six plateaus coincide, so that the overlap at the top of
    # the panel is quantified instead of merely looking like a single curve
    inner = np.abs(grids[REFERENCE][0]) <= 0.5
    stack = np.array([grids[step][1][inner] for step in steps])
    spread = float((stack.max(axis=0) - stack.min(axis=0)).max())
    axis.text(-0.66, 0.075,
              "on $|x| \\leq 0.5$ all six agree to "
              f"${spread * 1e4:.1f}\\times 10^{{-4}}$",
              fontsize=10, color="#6E6E6E", ha="left", va="center",
              style="italic")
    axis.text(-0.66, 0.017,
              r"$\Delta t = 0.4$ dashed: it coincides with $\Delta t = 0.3$",
              fontsize=10, color="#6E6E6E", ha="left", va="center",
              style="italic")
    print("plateau spread on |x|<=0.5:", spread)
    axis.set_title("(a)   final profiles at the common time $T = 19.2$",
                   fontsize=14, fontweight="bold", loc="left", pad=10)


def convergence_panel(axis, steps, rel, ceiling, order):
    series = (("inf", "#661100", r"$L^\infty$ error of $c$", "s"),
              ("l2", "#332288", r"$L^2$ error of $c$", "o"),
              ("front", "#117733", r"front position error", "D"))
    xlo, xhi = 0.0205, 0.58

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.add_patch(Rectangle((FILLED_FROM, 1e-4), xhi - FILLED_FROM, 100.0,
                             facecolor=BAND, edgecolor="none", zorder=0))
    axis.grid(True, which="major", linewidth=0.7, color="#DEDEDE", zorder=1)
    axis.set_axisbelow(False)

    for key, colour, label, marker in series:
        axis.plot([xlo, xhi], [ceiling[key]] * 2, ":", color=colour,
                  linewidth=1.8, alpha=0.75, zorder=3)
        xs = [s for s in steps if s != REFERENCE]
        ys = [rel[key][s] for s in xs]
        axis.plot(xs, ys, "-", color=colour, linewidth=3.0, zorder=4,
                  solid_capstyle="round")
        for s, y in zip(xs, ys):
            inside = s < FILLED_FROM
            axis.plot([s], [y], marker, color=colour, markersize=9,
                      markerfacecolor=colour if inside else "white",
                      markeredgewidth=2.2, zorder=5)

    axis.set_xlim(xlo, xhi)
    axis.set_ylim(0.035, 2.6)
    axis.set_xticks([0.025, 0.05, 0.1, 0.2, 0.4])
    axis.set_xticklabels(["0.025", "0.05", "0.1", "0.2", "0.4"])
    axis.set_yticks([0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
    axis.set_yticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2"])
    axis.tick_params(labelsize=13, which="major")
    axis.tick_params(which="minor", length=0)
    axis.set_xlabel(r"time step $\Delta t$", fontsize=15, fontweight="bold",
                    labelpad=1)
    axis.set_ylabel(r"relative error against the $\Delta t = 0.025$ run",
                    fontsize=13.5, fontweight="bold", labelpad=2)
    axis.set_title("(b)   convergence, and where it saturates",
                   fontsize=14, fontweight="bold", loc="left", pad=10)

    # direct labels, placed off the segment each one names
    places = {"inf": (0.0705, 1.14, "bottom"),
              "l2": (0.0705, 0.245, "bottom"),
              "front": (0.0705, 0.083, "top")}
    for key, colour, label, marker in series:
        x, y, va = places[key]
        axis.text(x, y, f"{label}    $p = {order[key]:.2f}$", fontsize=12,
                  color=colour, fontweight="bold", ha="center", va=va,
                  zorder=6)

    axis.text(np.sqrt(FILLED_FROM * xhi), 0.043,
              r"$\Delta t \geq 0.2$:" "\n"
              "no front anywhere,\n"
              "the interval is full\n"
              "(open markers)",
              fontsize=11, ha="center", va="bottom", color="#8A4A10",
              fontweight="bold", zorder=6, linespacing=1.4)
    axis.text(0.0225, 2.45,
              "dotted: the value each measure takes for\n"
              r"the trivial filled state $c \equiv 1$",
              fontsize=10.5, ha="left", va="top", color="#555555",
              zorder=6, linespacing=1.4)
    axis.annotate("reference run\n" r"error $\equiv 0$",
                  xy=(REFERENCE, 0.0365), xytext=(0.0222, 0.050),
                  fontsize=10.5, color="#555555", ha="left", va="bottom",
                  arrowprops=dict(arrowstyle="-|>", color="#999999",
                                  linewidth=1.4, shrinkA=3, shrinkB=0),
                  zorder=6, linespacing=1.4)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    return parser.parse_args()


def main():
    global RESULTS
    args = arguments()
    RESULTS = args.profiles.parent
    grids, summary = load()
    steps = sorted(grids)
    rel, fronts, minima, ceiling, order, front_ref = diagnostics(grids, steps)

    print("front positions", {k: (round(v, 7) if v else None)
                              for k, v in fronts.items()})
    print("minima", minima)
    for k in rel:
        print(k, {s: round(v, 6) for s, v in rel[k].items()},
              "ceiling", round(ceiling[k], 6), "order", round(order[k], 3))
    for step in steps:
        assert abs(summary[step]["max_error"]
                   - rel["inf"][step]) < 1e-12, step

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 12.5,
        "axes.linewidth": 2.0,
        "axes.edgecolor": INK,
        "xtick.color": INK, "ytick.color": INK, "text.color": INK,
        "xtick.major.width": 2.0, "ytick.major.width": 2.0,
        "xtick.major.size": 6.5, "ytick.major.size": 6.5,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })
    figure = plt.figure(figsize=(13.6, 6.4))
    grid = figure.add_gridspec(1, 2, width_ratios=(1.0, 1.0),
                               left=0.052, right=0.988, bottom=0.215,
                               top=0.930, wspace=0.20)
    left = figure.add_subplot(grid[0, 0])
    right = figure.add_subplot(grid[0, 1])
    for axis in (left, right):
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    profiles_panel(left, grids, steps, fronts, minima)
    convergence_panel(right, steps, rel, ceiling, order)

    figure.text(0.5, 0.088,
                "Fisher-Kolmogorov on $[-1,1]$, $\\alpha = 2$, $d = 0.0002$, "
                "200 P1 elements, backward Euler with Newton; all six runs "
                "stopped at $T = 19.2$, the largest time divisible by every "
                "step tested.",
                fontsize=10, color="#4A4A4A", ha="center", va="bottom")
    figure.text(0.5, 0.052,
                "In (b) each measure is divided by its own reference value, so "
                "all three are dimensionless and may share one axis; $p$ is "
                "the slope through $\\Delta t = 0.05$ and $0.1$ only.",
                fontsize=10, color="#4A4A4A", ha="center", va="bottom")
    figure.text(0.5, 0.016,
                "A front that has left the interval is censored at the domain "
                "end $x = 1$, so the three open diamonds are lower bounds and "
                "not measured front positions.",
                fontsize=10, color="#4A4A4A", ha="center", va="bottom")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=160, facecolor="white",
                   bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), facecolor="white",
                   bbox_inches="tight")
    print("saved", args.output)


if __name__ == "__main__":
    main()
