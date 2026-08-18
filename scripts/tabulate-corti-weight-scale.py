#!/usr/bin/env python3

"""Does the regional ranking of the deterministic model depend on the
undeclared scale of the connectivity weights?

Corti et al. do not state the scale of their edge weights, and benchmark 21
normalizes them as D_e = w_e / max(w), a choice made here. This script
tabulates, for runs of test_fisher_kolmogorov_corti83 at several multipliers
of that normalization, the regional ranking at the final time and its
number of pairwise agreements with table 3 of the reference, so that the
dependence of the reproduced ordering on the choice is on record.
"""

import argparse
import csv
import itertools
from pathlib import Path

GROUPS = ("frontal", "temporal", "parietal", "insular", "limbic",
          "occipital", "subcortical")
CORTI_RANK = {"frontal": 1, "limbic": 2, "temporal": 3, "insular": 4,
              "subcortical": 5, "parietal": 6, "occipital": 7}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("output/fisher_kolmogorov/corti83_scale"))
    parser.add_argument("--scales", nargs="+",
                        default=["35.3221", "10", "3", "1", "0.3", "0.1"])
    parser.add_argument("--maximum-weight", type=float, default=35.3221)
    parser.add_argument("--fiedler", type=float, default=0.772254)
    parser.add_argument("--mean-rate", type=float, default=0.1252)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    args = arguments()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["weight_scale", "rho", "damkohler",
                         "agreeing_pairs_of_21", "ranking_at_final_time"])
        for scale in args.scales:
            rows = list(csv.DictReader(
                open(args.root / f"s_{scale}" / "regional_averages.csv",
                     newline="")))
            values = {group: float(rows[-1][group]) for group in GROUPS}
            order = sorted(GROUPS, key=lambda group: -values[group])
            ours = {group: k + 1 for k, group in enumerate(order)}
            agree = sum(
                1 for a, b in itertools.combinations(GROUPS, 2)
                if (CORTI_RANK[a] < CORTI_RANK[b]) == (ours[a] < ours[b]))
            rho = float(scale) / args.maximum_weight
            damkohler = args.mean_rate / (rho * args.fiedler)
            writer.writerow([scale, f"{rho:.6g}", f"{damkohler:.4f}", agree,
                             " > ".join(order)])
            print(f"scale {scale:>8s}  rho {rho:8.4f}  Da {damkohler:6.2f}  "
                  f"agreement {agree:2d}/21  " + " > ".join(order))


if __name__ == "__main__":
    main()
