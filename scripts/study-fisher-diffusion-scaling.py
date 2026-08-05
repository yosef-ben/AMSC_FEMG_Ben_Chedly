#!/usr/bin/env python3

"""Diffusion-scaling study for the 83-region Fisher-Kolmogorov connectome.

Fornari et al. write the network model as dc/dt = -L c + alpha c (1-c) with L
the connectivity-weighted graph Laplacian. Their published lobe biomarkers are
clearly separated in time, whereas the literal Laplacian of the reconstructed
Budapest-83 graph homogenises the network long before the reaction becomes
visible. This script sweeps a uniform scaling rho of the connectivity weights,
measures the resulting lobe separation, and reports it against the Damkohler
number Da = alpha / (rho * lambda_2), lambda_2 being the Fiedler value of L.

It also records the largest rho at which the metric-graph FEM keeps the
concentration inside [0,1], which is the validity boundary of the consistent
mass P1 discretization used in benchmarks 19 and 21.
"""

import argparse
import csv
import json
import subprocess
from pathlib import Path

import numpy as np

DATA = Path("data/connectome/fornari83")
LOBES = ("temporal", "frontal", "parietal", "occipital")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build-release/test_fisher_kolmogorov_fornari83"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--dt", type=float, default=0.4)
    parser.add_argument("--final-time", type=float, default=120.0)
    parser.add_argument("--scalings", type=float, nargs="+",
                        default=[1.0, 0.5, 0.2, 0.1, 0.05, 0.02,
                                 0.01, 0.005, 0.002, 0.001])
    return parser.parse_args()


def graph_laplacian():
    edges = list(csv.DictReader(open(DATA / "edges.csv")))
    nodes = list(csv.DictReader(open(DATA / "nodes.csv")))
    n = len(nodes)
    adjacency = np.zeros((n, n))
    for row in edges:
        i, j = int(row["source"]), int(row["target"])
        w = float(row["connectivity_weight"])
        adjacency[i, j] += w
        adjacency[j, i] += w
    return np.diag(adjacency.sum(axis=1)) - adjacency, adjacency


def crossing_time(times, curve, level=50.0):
    """Linear interpolation of the first upward crossing of `level`."""
    for k in range(1, len(curve)):
        if curve[k - 1] < level <= curve[k]:
            span = curve[k] - curve[k - 1]
            weight = (level - curve[k - 1]) / span
            return times[k - 1] + weight * (times[k] - times[k - 1])
    return float("nan")


def run(executable, output_dir, scaling, alpha, dt, final_time, scheme):
    command = [str(executable), "1", str(dt), str(output_dir), str(alpha),
               str(final_time), str(scaling), scheme]
    return subprocess.run(command, text=True, capture_output=True)


def read_biomarkers(path):
    table = {name: [] for name in ("time", "global") + LOBES}
    with open(path) as stream:
        for row in csv.DictReader(stream):
            for name in table:
                table[name].append(float(row[name]))
    return {name: np.array(values) for name, values in table.items()}


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    laplacian, adjacency = graph_laplacian()
    spectrum = np.linalg.eigvalsh(laplacian)
    fiedler = spectrum[1]
    maximum_weight = adjacency.max()

    summary = {
        "fiedler_value": fiedler,
        "largest_eigenvalue": spectrum[-1],
        "mean_weighted_degree": float(np.diag(laplacian).mean()),
        "maximum_adjacency": float(maximum_weight),
        "benchmark_19_scaling": 1.0,
        "benchmark_21_scaling": float(1.0 / maximum_weight),
        "alpha": args.alpha,
        "dt": args.dt,
        "final_time": args.final_time,
    }
    print(f"Fiedler value lambda_2 = {fiedler:.6f}, "
          f"lambda_max = {spectrum[-1]:.4f}")
    print(f"benchmark 21 uses rho = 1/max(w) = {1.0 / maximum_weight:.6f}")

    rows = []
    curves = {}
    for scaling in args.scalings:
        work = args.output_dir / "runs" / f"rho_{scaling:g}"
        completed = run(args.executable, work, scaling, args.alpha,
                        args.dt, args.final_time, "nodal")
        if completed.returncode != 0:
            raise RuntimeError(f"nodal run failed at rho={scaling}:\n"
                               f"{completed.stderr}")
        table = read_biomarkers(work / "nodal_biomarkers.csv")
        curves[scaling] = table
        times = {name: crossing_time(table["time"], table[name])
                 for name in LOBES + ("global",)}
        lobe_times = [times[name] for name in LOBES]
        spread = max(lobe_times) - min(lobe_times)
        order = "<".join(sorted(LOBES, key=lambda name: times[name]))

        # Metric-graph FEM stability at one element per edge.
        status = {}
        for scheme in ("be", "cn"):
            result = run(args.executable, work / scheme, scaling, args.alpha,
                         args.dt, args.final_time, scheme)
            status[scheme] = "bounded" if result.returncode == 0 else "unbounded"

        row = {
            "diffusion_scaling": scaling,
            "damkohler": args.alpha / (scaling * fiedler),
            "t50_temporal": times["temporal"],
            "t50_frontal": times["frontal"],
            "t50_parietal": times["parietal"],
            "t50_occipital": times["occipital"],
            "t50_global": times["global"],
            "lobe_spread_years": spread,
            "activation_order": order,
            "fem_backward_euler": status["be"],
            "fem_corti_crank_nicolson": status["cn"],
        }
        rows.append(row)
        print(f"rho={scaling:<8g} Da={row['damkohler']:<9.3f} "
              f"spread={spread:8.4f} yr  FEM be/cn: "
              f"{status['be']}/{status['cn']}")

    with open(args.output_dir / "diffusion_scaling.csv", "w",
              newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    bounded = [row["diffusion_scaling"] for row in rows
               if row["fem_corti_crank_nicolson"] == "bounded"]
    summary["smallest_bounded_fem_scaling"] = min(bounded) if bounded else None
    with open(args.output_dir / "diffusion_scaling_summary.json", "w") as stream:
        json.dump(summary, stream, indent=2)

    np.savez(args.output_dir / "curves.npz",
             **{f"rho_{scaling:g}_{name}": table[name]
                for scaling, table in curves.items()
                for name in table})
    print(f"\nWritten to {args.output_dir}")


if __name__ == "__main__":
    main()
