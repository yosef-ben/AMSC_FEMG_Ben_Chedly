#!/usr/bin/env python3

"""Lobe separation at unit weight scale under every discretization we ran.

The comparison with Fornari et al. is made at the literal scale of the
connectivity weights, where the four lobe biomarkers coincide. This table
records that the coincidence is not an artefact of one discretization: it
collects, from the stored biomarker files, the 50-percent crossing spread of
the four lobes for the nodal model, for the metric-graph FEM with the fully
implicit and with the semi-implicit scheme, and for the FEM at eight elements
per connection over a fourfold-halved time step. The final range of each run
is stored next to it, so a run that ends outside [0,1] is visible as such
rather than dropped.
"""

import argparse
import csv
from pathlib import Path

LOBES = ("temporal", "frontal", "parietal", "occipital")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path,
                        default=Path("output/fisher_kolmogorov/fornari83"))
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
    network = crossing(times, [float(row["global"]) for row in rows])
    return (max(crossings) - min(crossings), network,
            float(rows[-1]["min"]), float(rows[-1]["max"]))


def main():
    args = arguments()
    base = args.base
    runs = [
        ("nodal", "backward_euler", 0.4, 1, base / "nodal_biomarkers.csv"),
        ("fem", "backward_euler", 0.4, 1, base / "fem_biomarkers.csv"),
        ("fem", "corti_semi_implicit", 0.4, 1,
         base / "scheme_cn_dt_0p4/fem_biomarkers.csv"),
        ("fem", "corti_semi_implicit", 0.4, 8,
         base / "scheme_cn_cells_8_dt_0p4/fem_biomarkers.csv"),
    ]
    for tag, step in (("0p8", 0.8), ("0p4", 0.4), ("0p2", 0.2),
                      ("0p1", 0.1), ("0p05", 0.05)):
        runs.append(("fem", "backward_euler", step, 8,
                     base / f"refinement/cells_8_dt_{tag}/fem_biomarkers.csv"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["model", "scheme", "dt", "cells_per_edge",
                         "lobe_spread_years", "t50_network_years",
                         "final_min", "final_max"])
        for model, scheme, step, cells, path in runs:
            spread, network, low, high = measure(path)
            writer.writerow([model, scheme, step, cells, f"{spread:.6g}",
                             f"{network:.6f}", f"{low:.6g}", f"{high:.6g}"])
            print(f"{model:5s} {scheme:20s} dt={step:<5g} cells={cells}  "
                  f"spread={spread:.3e} yr  t50={network:.4f}  "
                  f"final=[{low:.6g}, {high:.6g}]")


if __name__ == "__main__":
    main()
