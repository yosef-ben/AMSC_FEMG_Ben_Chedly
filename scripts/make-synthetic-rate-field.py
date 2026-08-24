#!/usr/bin/env python3

"""A synthetic per-vertex reaction field, for the control of benchmark 27.

The regional rates of Corti et al. are estimated for amyloid-beta, so using
them for the tau seeding is a transfer across proteins. This script writes a
rate field that shares nothing with them except the spread: it is linear in
the anterior-posterior coordinate of the vertices, fastest at the anterior
end, and its fastest vertex is 3.3 times its slowest, the ratio between the
frontal and the occipital coefficient of the reference. The executable
rescales whatever it reads to the requested vertex mean, so only the shape of
the field enters. The file is written in the format the solver reads, the
same one test_fisher_kolmogorov_corti83 produces.
"""

import argparse
import csv
from pathlib import Path

NODES = Path("data/connectome/fornari83/nodes.csv")
# Frontal over occipital coefficient of table 1 of Corti et al.
RATIO = 0.1801 / 0.0545


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path, default=NODES)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    args = arguments()
    with open(args.nodes, newline="") as stream:
        nodes = list(csv.DictReader(stream))
    anterior = [float(node["y"]) for node in nodes]
    low, high = min(anterior), max(anterior)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["node_id", "name", "region", "alpha"])
        for node in nodes:
            fraction = (float(node["y"]) - low) / (high - low)
            rate = 1.0 + (RATIO - 1.0) * fraction
            writer.writerow([node["node_id"], node["name"], "synthetic",
                             f"{rate:.10f}"])
    print(f"Written {args.output}: linear in y, ratio {RATIO:.4f} between "
          f"the anterior and the posterior end")


if __name__ == "__main__":
    main()
