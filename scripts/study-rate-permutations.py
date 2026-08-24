#!/usr/bin/env python3

"""What a regional reaction field must satisfy to reorder the lobes.

The seven coefficients this chapter borrows are estimated by Corti et al. for
amyloid-beta, so whether they are compatible with tau is not established. This
study removes the question from the conclusion: the seven values are assigned
to the seven anatomical groups in all 5040 possible ways, each assignment is
rescaled to the same vertex mean, and the resulting order of the four lobe
biomarkers is recorded. What the clinical sequence needs is then read off the
table rather than assumed.

The model is the one of benchmark 27, the metric-graph FEM with the lumped
mass at one element per connection, written as the ordinary differential
system it reduces to, with the entorhinal seed and the transport scaling
rho = 0.005. Heun's method with a small step is used instead of the backward
Euler of the solver, so the crossing times are the small-step limit and are
some two years earlier than the stored runs at dt = 0.4; only the order of the
four lobes is used here, and it is the same in both.
"""

import argparse
import csv
import itertools
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lobe_scale import classify

NODES = Path("data/connectome/fornari83/nodes.csv")
EDGES = Path("data/connectome/fornari83/edges.csv")
REGIONS = Path("benchmarks/21_fisher_kolmogorov_corti83/results"
               "/reaction_coefficients.csv")
GROUPS = ("frontal", "temporal", "parietal", "insular", "limbic",
          "occipital", "subcortical")
# Table 1 of Corti et al., in the order of GROUPS.
VALUES = np.array([0.1801, 0.1421, 0.0627, 0.1005, 0.1351, 0.0545, 0.1147])
LOBES = ("temporal", "frontal", "parietal", "occipital")
CLINICAL = ("temporal", "frontal", "parietal", "occipital")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho", type=float, default=0.005)
    parser.add_argument("--mean-rate", type=float, default=0.5)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--final-time", type=float, default=90.0)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_csv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    args = arguments()
    nodes = {int(r["node_id"]): r["name"] for r in read_csv(NODES)}
    size = 1 + max(nodes)
    lobe = {k: classify(name) for k, name in nodes.items()}
    group = {int(r["node_id"]): r["region"] for r in read_csv(REGIONS)}

    laplacian = np.zeros((size, size))
    lumped = np.zeros(size)
    for edge in read_csv(EDGES):
        i, j = int(edge["source"]), int(edge["target"])
        weight = float(edge["connectivity_weight"])
        laplacian[i, i] += weight
        laplacian[j, j] += weight
        laplacian[i, j] -= weight
        laplacian[j, i] -= weight
        lumped[i] += 0.5
        lumped[j] += 0.5
    transport = (args.rho / lumped)[:, None] * laplacian
    seed = [k for k, name in nodes.items() if "entorhinal" in name.lower()]
    in_group = {g: np.array([k for k in range(size) if group[k] == g])
                for g in GROUPS}
    in_lobe = {l: np.array([k for k in range(size) if lobe[k] == l])
               for l in LOBES}

    permutations = list(itertools.permutations(range(len(GROUPS))))
    count = len(permutations)
    rates = np.zeros((count, size))
    for p, permutation in enumerate(permutations):
        for position, g in enumerate(GROUPS):
            rates[p, in_group[g]] = VALUES[permutation[position]]
    rates *= args.mean_rate * size / rates.sum(axis=1)[:, None]

    solution = np.zeros((count, size))
    solution[:, seed] = 0.1
    crossing = {l: np.full(count, np.nan) for l in LOBES}
    previous = {l: solution[:, in_lobe[l]].mean(axis=1) for l in LOBES}
    time = 0.0

    def rate_of_change(state):
        return -(state @ transport.T) + rates * state * (1.0 - state)

    while time < args.final_time:
        first = rate_of_change(solution)
        second = rate_of_change(solution + args.dt * first)
        solution = solution + 0.5 * args.dt * (first + second)
        time += args.dt
        for l in LOBES:
            current = solution[:, in_lobe[l]].mean(axis=1)
            hit = (np.isnan(crossing[l]) & (previous[l] < 0.5)
                   & (current >= 0.5))
            crossing[l][hit] = time - args.dt + args.dt * (
                0.5 - previous[l][hit]) / (current[hit] - previous[l][hit])
            previous[l] = current

    times = np.vstack([crossing[l] for l in LOBES]).T
    complete = ~np.isnan(times).any(axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    clinical = 0
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["assignment"] + [f"alpha_{g}" for g in GROUPS]
                        + ["frontal_over_occipital", "order",
                           "clinical_order"])
        for p in range(count):
            if not complete[p]:
                order = "incomplete"
                is_clinical = 0
            else:
                order = " ".join(LOBES[i] for i in np.argsort(times[p]))
                is_clinical = int(order == " ".join(CLINICAL))
            clinical += is_clinical
            ratio = (rates[p, in_group["frontal"]][0]
                     / rates[p, in_group["occipital"]][0])
            writer.writerow([p] + [f"{rates[p, in_group[g]][0]:.6f}"
                                   for g in GROUPS]
                            + [f"{ratio:.6f}", order, is_clinical])
    print(f"Written {args.output}: {count} assignments, {clinical} give the "
          f"clinical order ({100.0 * clinical / count:.1f} percent)")


if __name__ == "__main__":
    main()
