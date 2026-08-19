#!/usr/bin/env python3

"""The Damkohler number of the nodal model against that of the metric-graph FEM.

Both models share the connectivity-weighted Laplacian L of the 83-region
graph, which the finite element diffusion matrix reproduces exactly at one
element per connection with unit lengths. They differ in the mass matrix: the
nodal model has the identity, the finite element model gives every vertex the
mass of its incident half-connections, M_ii = deg_i / 2 with the lumped mass
(h/2 per node and cell) and the consistent blocks h/3, h/6 otherwise. The rate
of the slowest transport mode is therefore the Fiedler value of L for the
nodal model and the smallest non-zero generalized eigenvalue of (L, M) for the
finite element one, and the Damkohler number alpha / (rho lambda_2) of the two
models differs by their ratio at the same rho. This script computes the three
values from the stored edge list and writes them next to the sweep results.
"""

import argparse
import csv
from pathlib import Path

import numpy as np


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", type=Path,
                        default=Path("data/connectome/fornari83/edges.csv"))
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--scalings", type=float, nargs="+",
                        default=[1.0, 0.05, 0.005])
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/23_fisher_kolmogorov_diffusion_scaling/results"
        "/mass_spectrum.csv"))
    return parser.parse_args()


def matrices(edges_path):
    with open(edges_path, newline="") as stream:
        edges = list(csv.DictReader(stream))
    n = 1 + max(max(int(e["source"]), int(e["target"])) for e in edges)
    laplacian = np.zeros((n, n))
    lumped = np.zeros(n)
    consistent = np.zeros((n, n))
    for edge in edges:
        i, j = int(edge["source"]), int(edge["target"])
        w = float(edge["connectivity_weight"])
        laplacian[i, i] += w
        laplacian[j, j] += w
        laplacian[i, j] -= w
        laplacian[j, i] -= w
        # One P1 element of unit length per connection.
        lumped[i] += 0.5
        lumped[j] += 0.5
        consistent[i, i] += 1.0 / 3.0
        consistent[j, j] += 1.0 / 3.0
        consistent[i, j] += 1.0 / 6.0
        consistent[j, i] += 1.0 / 6.0
    return laplacian, lumped, consistent


def fiedler(laplacian, mass=None):
    """Smallest non-zero eigenvalue of L v = lambda M v, M symmetric positive
    definite (the identity when omitted)."""
    if mass is None:
        reduced = laplacian
    else:
        factor = np.linalg.cholesky(mass)
        inverse = np.linalg.inv(factor)
        reduced = inverse @ laplacian @ inverse.T
    values = np.sort(np.linalg.eigvalsh(reduced))
    return float(values[1])


def main():
    args = arguments()
    laplacian, lumped, consistent = matrices(args.edges)
    rows = [
        ("identity", fiedler(laplacian), 1.0, 1.0, 1.0),
        ("lumped", fiedler(laplacian, np.diag(lumped)),
         float(lumped.mean()), float(lumped.min()), float(lumped.max())),
        ("consistent", fiedler(laplacian, consistent),
         float(np.diag(consistent).mean()), float(np.diag(consistent).min()),
         float(np.diag(consistent).max())),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mass_matrix", "lambda_2", "mean_diagonal",
                         "min_diagonal", "max_diagonal"]
                        + [f"damkohler_rho_{s:g}" for s in args.scalings])
        for name, value, mean, low, high in rows:
            writer.writerow([name, f"{value:.6g}", f"{mean:.6g}",
                             f"{low:.6g}", f"{high:.6g}"]
                            + [f"{args.alpha / (s * value):.6g}"
                               for s in args.scalings])
    reference = rows[0][1]
    for name, value, mean, low, high in rows:
        print(f"{name:11s} lambda_2 = {value:.5f}  ratio to identity = "
              f"{reference / value:5.2f}  diagonal mean {mean:5.2f} "
              f"[{low:g}, {high:g}]  Da at rho = "
              + ", ".join(f"{s:g}: {args.alpha / (s * value):.4g}"
                          for s in args.scalings))
    print(f"Written {args.output}")


if __name__ == "__main__":
    main()
