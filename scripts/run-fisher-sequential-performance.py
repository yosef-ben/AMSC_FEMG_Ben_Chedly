#!/usr/bin/env python3

"""Measure the sequential Corti-83 baseline without visualization I/O."""

import argparse
import csv
import statistics
import subprocess
import os
from pathlib import Path


FIELDS = [
    "n_vertices", "n_edges", "n_dofs", "matrix_nnz",
    "init_seconds", "assembly_seconds", "solve_seconds", "total_seconds",
]


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build/test_fisher_kolmogorov_corti83"))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def run(executable, cells):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    completed = subprocess.run(
        [str(executable), str(graph), "--performance"],
        check=True, text=True, capture_output=True, env=environment,
    )
    line = next((line for line in completed.stdout.splitlines()
                 if line.startswith("PERFORMANCE,")), None)
    if line is None:
        raise RuntimeError(f"Missing timing record:\n{completed.stdout}")
    values = line.split(",")[1:]
    record = {name: float(value) for name, value in zip(FIELDS, values)}
    for name in ("n_vertices", "n_edges", "n_dofs", "matrix_nnz"):
        record[name] = int(record[name])
    record["cells_per_edge"] = cells
    return record


def main():
    args = arguments()
    if args.repeats < 1 or args.warmups < 0:
        raise ValueError("Invalid repeat count.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = []
    summary = []
    for cells in (1, 2, 4, 8):
        for _ in range(args.warmups):
            run(args.executable, cells)
        samples = []
        for repeat in range(args.repeats):
            record = run(args.executable, cells)
            record["repeat"] = repeat + 1
            samples.append(record)
            raw.append(record)
        row = {key: samples[0][key] for key in
               ("cells_per_edge", "n_vertices", "n_edges",
                "n_dofs", "matrix_nnz")}
        for key in ("init_seconds", "assembly_seconds",
                    "solve_seconds", "total_seconds"):
            values = [sample[key] for sample in samples]
            row[key] = statistics.median(values)
            row[key.replace("_seconds", "_min_seconds")] = min(values)
            row[key.replace("_seconds", "_max_seconds")] = max(values)
        row["seconds_per_step"] = row["solve_seconds"] / 100.0
        summary.append(row)
        print(f"cells={cells}: DoFs={row['n_dofs']}, "
              f"solve={row['solve_seconds']:.6f} s")

    raw_fields = ["cells_per_edge", "repeat", *FIELDS]
    with (args.output_dir / "raw_timings.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=raw_fields)
        writer.writeheader()
        writer.writerows(raw)

    summary_fields = list(summary[0])
    with (args.output_dir / "sequential_performance.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
