#!/usr/bin/env python3

"""Lobe separation of the consistent-mass FEM within its validity boundary.

The consistent mass matrix keeps the solution bounded only up to the scaling
rho = 0.05 (benchmark 23), so its sweep cannot be tabulated beyond that
point. This script reads the stored runs of the fully implicit scheme at one
element per connection for the bounded scalings and records, for each, the
50-percent crossing spread of the four lobe biomarkers, the transient
extremes and the network crossing, in the format of fem_lumped_sweep.csv.
"""

import argparse
import csv
from pathlib import Path

FIEDLER = 0.772254
LOBES = ("temporal", "frontal", "parietal", "occipital")


def crossing(times, values, level=50.0):
    for k in range(1, len(values)):
        if values[k - 1] < level <= values[k]:
            span = values[k] - values[k - 1]
            return times[k - 1] + (level - values[k - 1]) / span * (
                times[k] - times[k - 1])
    return float("nan")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("output/fisher_kolmogorov/stabilization"))
    parser.add_argument("--scalings", nargs="+",
                        default=["1.0", "0.5", "0.2", "0.1", "0.05"])
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/23_fisher_kolmogorov_diffusion_scaling/results"
        "/fem_consistent_bounded.csv"))
    return parser.parse_args()


def main():
    args = arguments()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["diffusion_scaling", "damkohler", "cells_per_edge",
                         "scheme", "transient_min", "transient_max",
                         "lobe_spread_years", "t50_network_years"])
        for scaling in args.scalings:
            path = args.root / f"rho_{scaling}_be" / "fem_biomarkers.csv"
            with open(path, newline="") as run:
                rows = list(csv.DictReader(run))
            times = [float(r["time"]) for r in rows]
            crossings = [crossing(times, [float(r[l]) for r in rows])
                         for l in LOBES]
            writer.writerow([
                scaling, f"{args.alpha / (float(scaling) * FIEDLER):.4f}", 1,
                "be",
                f"{min(float(r['min']) for r in rows):.6g}",
                f"{max(float(r['max']) for r in rows):.6g}",
                f"{max(crossings) - min(crossings):.6g}",
                f"{crossing(times, [float(r['global']) for r in rows]):.6f}"])
            print(f"  rho = {scaling:>5s}  spread = "
                  f"{max(crossings) - min(crossings):.3f} yr")
    print(f"Written {args.output}")


if __name__ == "__main__":
    main()
