#!/usr/bin/env python3

"""Instability of the consistent-mass FEM and its cure by mass lumping.

Two rows of three panels, all at one element per connection with the fully
implicit scheme and dt = 0.4 years, at three transport scalings. Top: the
consistent mass matrix of the first part of the report; bottom: the row-sum
lumped mass with the vertex-rule reaction. Each panel draws the four lobe
biomarkers in the colours of Fornari et al. and, shaded, the envelope of the
83 vertex concentrations, from the smallest to the largest. Where the
consistent-mass run leaves the physical range the envelope crosses the dashed
lines at 0 and 100 percent; at the weakest transport the Newton iteration
fails and the curves stop, which is drawn as it happened. Nothing is clipped
or rescaled; the axis is the same in every panel.
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

LOBE_STYLE = {
    "temporal": ("#2CA02C", "-"),
    "frontal": ("#D62728", "--"),
    "parietal": ("#FF7F0E", "-."),
    "occipital": ("#1F77B4", ":"),
}
ROWS = (("be", "consistent mass"), ("be_lumped", "lumped mass"))


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("output/fisher_kolmogorov/stabilization"))
    parser.add_argument("--scalings", nargs="+", default=["1.0", "0.05", "0.005"])
    parser.add_argument("--fiedler", type=float, default=0.772254)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read(path):
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {key: [float(row[key]) for row in rows] for key in rows[0]}


def main():
    args = arguments()
    figure_style.apply()
    figure, axes = plt.subplots(2, 3, figsize=(9.8, 6.0), sharex=True,
                                sharey=True)
    summary = []
    for row, (scheme, title) in enumerate(ROWS):
        for column, scaling in enumerate(args.scalings):
            axis = axes[row, column]
            data = read(args.root / f"rho_{scaling}_{scheme}"
                        / "fem_biomarkers.csv")
            time = data["time"]
            low = [100.0 * v for v in data["min"]]
            high = [100.0 * v for v in data["max"]]
            axis.fill_between(time, low, high, color="0.75", alpha=0.45,
                              linewidth=0, zorder=1)
            axis.axhline(0.0, color="0.5", linewidth=0.9,
                         linestyle=(0, (4, 3)), zorder=2)
            axis.axhline(100.0, color="0.5", linewidth=0.9,
                         linestyle=(0, (4, 3)), zorder=2)
            for lobe, (colour, style) in LOBE_STYLE.items():
                axis.plot(time, data[lobe], color=colour, linestyle=style,
                          linewidth=1.7, zorder=3)
            damkohler = args.alpha / (float(scaling) * args.fiedler)
            if row == 0:
                # Lower-right corner, empty in every panel of the top row.
                axis.text(0.62, 0.30,
                          rf"$\rho = {float(scaling):g}$" "\n"
                          rf"Da $= {damkohler:.1f}$",
                          transform=axis.transAxes, fontsize=9,
                          fontweight="bold", color="0.35", va="top")
            axis.set_xlim(0, 40)
            axis.set_ylim(-60, 115)
            axis.set_xticks([0, 10, 20, 30, 40])
            axis.set_yticks([-50, 0, 50, 100])
            letter = chr(97 + 3 * row + column)
            axis.text(0.0, 1.04, f"({letter})", transform=axis.transAxes,
                      fontsize=10.5, fontweight="bold", style="italic",
                      va="bottom")
            if column == 0:
                axis.set_ylabel(f"{title}\nbiomarker [%]", labelpad=2)
            if row == 1:
                figure_style.xname(axis, "t [yr]", y=-0.14, fontsize=9)
            summary.append((scheme, scaling, min(data["min"]),
                            max(data["max"]), time[-1]))
    for slot, (lobe, (colour, _)) in enumerate(LOBE_STYLE.items()):
        figure_style.label_series(axes[1, 0], 24.0, 44.0 - 13.0 * slot,
                                  lobe, colour, fontsize=8.5)
    figure.tight_layout(w_pad=1.6, h_pad=2.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)
    figure.savefig(args.output.with_suffix(".pdf"))
    print(f"Written {args.output} and its PDF")
    # The numbers the caption and the text quote, next to the figure.
    with open(args.output.with_name("stabilization_summary.csv"), "w",
              newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scheme", "diffusion_scaling", "damkohler",
                         "min_concentration", "max_concentration",
                         "last_stored_time"])
        for scheme, scaling, low, high, last in summary:
            writer.writerow([scheme, scaling,
                             f"{args.alpha / (float(scaling) * args.fiedler):.4f}",
                             f"{low:.6g}", f"{high:.6g}", f"{last:g}"])
    for scheme, scaling, low, high, last in summary:
        print(f"  {scheme:10s} rho={scaling:>6s}  min c = {low:9.4g}  "
              f"max c = {high:8.5f}  last stored time = {last:g}")


if __name__ == "__main__":
    main()
