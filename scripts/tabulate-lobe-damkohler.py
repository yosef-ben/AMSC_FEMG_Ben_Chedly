#!/usr/bin/env python3

"""Global and lobe-scale transport rates of the three models, as a record.

For the nodal model and for the metric-graph FEM with consistent and lumped
mass, the smallest non-zero eigenvalue of the Laplacian relative to the mass
matrix (the global rate behind the Damkohler number of the report), the
relaxation rates of the patterns constant on the four lobes and the
remaining regions (the compressed graph of the report, whose slowest rate is
the lobe-scale eigenvalue), the lobe-scale Damkohler number at the three
scalings of the report figure and the onset ratio between the lobe rate and
the nodal global rate. A second file records the largest components of the
eigenvector of the nodal global rate, which show that the slowest pattern of
the network is not a contrast between lobes.
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lobe_scale import LobeGraph, MODELS, damkohler_lobe


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--scalings", type=float, nargs="+",
                        default=[1.0, 0.05, 0.005])
    parser.add_argument("--output", type=Path, default=Path(
        "benchmarks/23_fisher_kolmogorov_diffusion_scaling/results"
        "/lobe_damkohler.csv"))
    return parser.parse_args()


def main():
    args = arguments()
    graph = LobeGraph()
    nodal_global = graph.global_rate("nodal")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["model", "global_rate", "lobe_rate", "lobe_rates",
                         "temporal_occipital_rate", "onset_ratio"]
                        + [f"damkohler_lobe_rho_{s:g}" for s in args.scalings])
        for model in MODELS:
            rates = graph.lobe_rates(model)
            writer.writerow([
                model, f"{graph.global_rate(model):.6g}", f"{rates[0]:.6g}",
                " ".join(f"{r:.6g}" for r in rates),
                f"{graph.contrast_rate('temporal', 'occipital', model):.6g}",
                f"{rates[0] / nodal_global:.6g}"]
                + [f"{damkohler_lobe(s, args.alpha, rates[0]):.6g}"
                   for s in args.scalings])
    support = args.output.with_name("fiedler_support.csv")
    with open(support, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["rank", "region", "component"])
        for rank, (name, value) in enumerate(graph.fiedler_support(), 1):
            writer.writerow([rank, name, f"{value:.6g}"])
    print(f"Written {args.output} and {support}")
    for model in MODELS:
        rates = graph.lobe_rates(model)
        print(f"  {model:10s} global {graph.global_rate(model):.4f}  lobe "
              f"{rates[0]:.3f}  onset ratio {rates[0] / nodal_global:.2f}  "
              "Da_lobe at rho = "
              + ", ".join(f"{s:g}: {damkohler_lobe(s, args.alpha, rates[0]):.3g}"
                          for s in args.scalings))


if __name__ == "__main__":
    main()
