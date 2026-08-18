#!/usr/bin/env python3

"""Boundedness of the mass-lumped finite element model over the sweep.

The consistent-mass P1 discretization leaves [0,1] beyond Da of about 13
(benchmark 23). Mass lumping, the row-sum diagonal mass with the reaction
evaluated by the vertex rule, is the standard remedy and is available in
``fisher_kolmogorov_problem`` through ``set_mass_lumping``. This script
tabulates, from the stored biomarker files of the lumped runs, the transient
extremes, the lobe spread and the network crossing at every scaling of the
sweep and for both time schemes, plus the four-element case at rho = 0.005
where the consistent mass diverged, so that the statement "the lumped
formulation stays within [0,1] over the whole sweep" rests on a stored table.
"""

import argparse
import csv
from pathlib import Path

LOBES = ("temporal", "frontal", "parietal", "occipital")
SCALINGS = ("1.0", "0.5", "0.2", "0.1", "0.05", "0.04", "0.03", "0.025",
            "0.02", "0.01", "0.005", "0.002", "0.001")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("output/fisher_kolmogorov/lumped"))
    parser.add_argument("--fiedler", type=float, default=0.772254)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def crossing(times, values, level=50.0):
    for k in range(1, len(values)):
        if values[k - 1] < level <= values[k]:
            span = values[k] - values[k - 1]
            return times[k - 1] + (level - values[k - 1]) / span * (
                times[k] - times[k - 1])
    return float("nan")


def measure(path):
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    times = [float(row["time"]) for row in rows]
    crossings = [crossing(times, [float(row[lobe]) for row in rows])
                 for lobe in LOBES]
    return (min(float(row["min"]) for row in rows),
            max(float(row["max"]) for row in rows),
            max(crossings) - min(crossings),
            crossing(times, [float(row["global"]) for row in rows]))


def main():
    args = arguments()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    runs = [(rho, 1, scheme) for rho in SCALINGS
            for scheme in ("be_lumped", "cn_lumped")]
    runs += [("0.005", 4, "be_lumped"), ("0.005", 4, "cn_lumped")]
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["diffusion_scaling", "damkohler", "cells_per_edge",
                         "scheme", "transient_min", "transient_max",
                         "lobe_spread_years", "t50_network_years"])
        for rho, cells, scheme in runs:
            folder = (args.root / f"rho_{rho}_{scheme}" if cells == 1
                      else args.root / f"rho_{rho}_cells{cells}_{scheme}")
            low, high, spread, network = measure(folder / "fem_biomarkers.csv")
            damkohler = args.alpha / (float(rho) * args.fiedler)
            writer.writerow([rho, f"{damkohler:.4f}", cells, scheme,
                             f"{low:.6g}", f"{high:.6g}", f"{spread:.6g}",
                             f"{network:.6f}"])
            print(f"rho={rho:>6s} Da={damkohler:7.2f} cells={cells} "
                  f"{scheme:9s} range=[{low:.2e}, {high:.6f}] "
                  f"spread={spread:.3e} t50={network:.4f}")


if __name__ == "__main__":
    main()
