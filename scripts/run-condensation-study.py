#!/usr/bin/env python3

"""Run the static-condensation study on the Corti-83 reference problem.

One process per mesh: the driver advances the semi-implicit scheme with the
per-edge static condensation and with the full-system SimplicialLDLT on the
natural ordering (the sequential reference of the ordering study), reports
the phase medians of both and the final-state difference between the two.
The validation records compare the condensed trajectory against the solve()
of the production class at every step. Threads pinned to one.
"""

import argparse
import csv
import os
import subprocess
from pathlib import Path

FIELDS = [
    "n_dofs", "steps",
    "cond_local_median_s", "cond_interface_median_s", "cond_back_median_s",
    "cond_step_median_s", "cond_loop_total_s",
    "full_rebuild_median_s", "full_factor_median_s", "full_solve_median_s",
    "full_step_median_s", "full_loop_total_s",
    "max_diff_condensed_vs_full",
]
CELLS = (8, 16, 32, 64, 128)
VALIDATION_CELLS = (1, 2, 4, 8)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build-release/test_condensation"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=100)
    return parser.parse_args()


def run(executable, cells, argument, max_steps):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    completed = subprocess.run(
        [str(executable), str(graph), argument, str(max_steps)],
        check=True, text=True, capture_output=True, env=environment,
    )
    return completed.stdout


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.output_dir / "validation.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", "n_dofs", "steps",
                         "max_diff_condensed_vs_class"])
        for cells in VALIDATION_CELLS:
            stdout = run(args.executable, cells, "--validate",
                         args.max_steps)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("VALIDATION,"))
            writer.writerow([cells, *line.split(",")[1:]])
            print(f"validate cells={cells}: max diff "
                  f"{line.split(',')[-1]}")

    with open(args.output_dir / "condensation_study.csv", "w",
              newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", *FIELDS])
        for cells in CELLS:
            stdout = run(args.executable, cells, "benchmark",
                         args.max_steps)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("CONDENSE,"))
            values = line.split(",")[1:]
            writer.writerow([cells, *values])
            record = dict(zip(FIELDS, values))
            print(f"cells={cells:>4} n={record['n_dofs']:>7} "
                  f"condensed {float(record['cond_step_median_s'])*1e3:7.2f} ms"
                  f" (local {float(record['cond_local_median_s'])*1e3:6.2f}"
                  f", iface {float(record['cond_interface_median_s'])*1e3:5.2f}"
                  f", back {float(record['cond_back_median_s'])*1e3:5.2f})"
                  f"  full {float(record['full_step_median_s'])*1e3:7.2f} ms"
                  f"  diff {record['max_diff_condensed_vs_full']}")


if __name__ == "__main__":
    main()
