#!/usr/bin/env python3

"""Run the sequential ordering study on the Corti-83 reference problem.

One process per (mesh, variant): the driver replays the production
semi-implicit loop with the requested solver-ordering variant and prints one
CSV record with the structural metrics (bandwidth, profile, factor nonzeros)
and the per-step phase times (matrix rebuild, factorization, numeric-only
refactorization, triangular solves). The validation records prove that the
replayed lu_colamd loop reproduces the solve() of the library class.

Threads are pinned to one everywhere: this is the sequential baseline.
"""

import argparse
import csv
import os
import subprocess
from pathlib import Path

FIELDS = [
    "variant", "n_dofs", "k_nnz", "bandwidth", "profile", "factor_nnz",
    "steps_timed", "failed", "rebuild_median_s", "factor_median_s",
    "refactor_median_s", "solve_median_s", "step_median_s", "loop_total_s",
    "max_diff_vs_colamd",
]
VARIANTS = ("lu_colamd", "lu_natural", "lu_amd",
            "ldlt_amd", "ldlt_natural", "lu_rcm")
CELLS = (1, 2, 4, 8, 16, 32, 64)
VALIDATION_CELLS = (1, 2, 4, 8)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build-release/test_ordering_study"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=100)
    parser.add_argument("--step-time-limit", type=float, default=5.0)
    return parser.parse_args()


def run(executable, cells, argument, max_steps, limit):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    completed = subprocess.run(
        [str(executable), str(graph), argument,
         str(max_steps), str(limit)],
        check=True, text=True, capture_output=True, env=environment,
    )
    return completed.stdout


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.output_dir / "validation.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", "n_dofs", "steps",
                         "max_diff_replayed_vs_class"])
        for cells in VALIDATION_CELLS:
            stdout = run(args.executable, cells, "--validate",
                         args.max_steps, args.step_time_limit)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("VALIDATION,"))
            values = line.split(",")[1:]
            writer.writerow([cells, *values])
            print(f"validate cells={cells}: max diff {values[-1]}")

    with open(args.output_dir / "ordering_study.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", *FIELDS])
        for cells in CELLS:
            for variant in VARIANTS:
                stdout = run(args.executable, cells, variant,
                             args.max_steps, args.step_time_limit)
                line = next(line for line in stdout.splitlines()
                            if line.startswith("ORDERING,"))
                values = line.split(",")[1:]
                writer.writerow([cells, *values])
                record = dict(zip(FIELDS, values))
                print(f"cells={cells:>3} {variant:<12} "
                      f"factor {float(record['factor_median_s'])*1e3:8.2f} ms"
                      f"  solve {float(record['solve_median_s'])*1e3:7.2f} ms"
                      f"  fill {record['factor_nnz']:>9}"
                      f"  steps {record['steps_timed']:>3}"
                      f"  diff {record['max_diff_vs_colamd']}")


if __name__ == "__main__":
    main()
