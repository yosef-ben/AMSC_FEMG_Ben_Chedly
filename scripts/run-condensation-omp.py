#!/usr/bin/env python3

"""Strong scaling of the OpenMP condensation on the 128-cell connectome.

One process per thread count, threads pinned (OMP_PLACES=cores,
OMP_PROC_BIND=close). Every benchmark process runs the sequential condensed
loop and the OpenMP loop back to back, so each record carries its own
in-process sequential reference; the speedup column uses the sequential
step median of the same row. The validation records advance the two
engines in lockstep and compare every step.
"""

import argparse
import csv
import os
import subprocess
from pathlib import Path

RAW_FIELDS = [
    "n_dofs", "threads", "steps",
    "seq_local_median_s", "seq_interface_median_s", "seq_back_median_s",
    "seq_step_median_s", "seq_loop_total_s",
    "omp_local_median_s", "omp_reduce_median_s", "omp_interface_median_s",
    "omp_back_median_s", "omp_step_median_s", "omp_loop_total_s",
    "max_diff_omp_vs_seq",
]
THREADS = (1, 2, 4, 8)
VALIDATIONS = ((8, 2), (8, 4), (8, 8), (128, 8))
SCALING_CELLS = 128
REPEATS = 3


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build-release/test_condensation_omp"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=200)
    return parser.parse_args()


def run(executable, cells, mode, threads, max_steps):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    environment = os.environ.copy()
    environment.update({
        "OMP_PLACES": "cores",
        "OMP_PROC_BIND": "close",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    completed = subprocess.run(
        [str(executable), str(graph), mode, str(threads), str(max_steps)],
        check=True, text=True, capture_output=True, env=environment,
    )
    return completed.stdout


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.output_dir / "omp_validation.csv", "w",
              newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", "n_dofs", "threads", "steps",
                         "max_diff_omp_vs_seq"])
        for cells, threads in VALIDATIONS:
            stdout = run(args.executable, cells, "--validate", threads,
                         args.max_steps)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("VALIDATION,"))
            writer.writerow([cells, *line.split(",")[1:]])
            print(f"validate cells={cells} threads={threads}: "
                  f"max diff {line.split(',')[-1]}")

    # Repeated processes absorb the turbo and scheduling noise of the
    # laptop CPU: every field is the median across the repeats and the
    # sequential reference T1 is pooled across every process.
    import statistics
    aggregated = {}
    for threads in THREADS:
        repeats = []
        for _ in range(REPEATS):
            stdout = run(args.executable, SCALING_CELLS, "benchmark",
                         threads, args.max_steps)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("OMP,"))
            repeats.append(dict(zip(RAW_FIELDS, line.split(",")[1:])))
        record = {}
        for name in RAW_FIELDS:
            if name in ("n_dofs", "threads", "steps"):
                record[name] = repeats[0][name]
            elif name == "max_diff_omp_vs_seq":
                record[name] = max(r[name] for r in repeats)
            else:
                record[name] = statistics.median(
                    float(r[name]) for r in repeats)
        aggregated[threads] = record
    t1 = statistics.median(
        [float(r["seq_step_median_s"]) for r in aggregated.values()])

    with open(args.output_dir / "omp_scaling.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", *RAW_FIELDS,
                         "t1_pooled_s", "speedup", "efficiency"])
        for threads in THREADS:
            record = aggregated[threads]
            speedup = t1 / float(record["omp_step_median_s"])
            efficiency = speedup / threads
            writer.writerow([SCALING_CELLS,
                             *[record[name] for name in RAW_FIELDS],
                             f"{t1:.6f}",
                             f"{speedup:.4f}", f"{efficiency:.4f}"])
            print(f"threads={threads}: "
                  f"seq {float(record['seq_step_median_s'])*1e3:6.2f} ms  "
                  f"omp {float(record['omp_step_median_s'])*1e3:6.2f} ms "
                  f"(local {float(record['omp_local_median_s'])*1e3:6.2f}"
                  f", reduce {float(record['omp_reduce_median_s'])*1e3:5.2f}"
                  f", iface {float(record['omp_interface_median_s'])*1e3:5.2f}"
                  f", back {float(record['omp_back_median_s'])*1e3:5.2f})  "
                  f"S={speedup:.2f} E={efficiency:.2f}  "
                  f"diff {record['max_diff_omp_vs_seq']}")
    print(f"pooled sequential T1: {t1*1e3:.2f} ms per step")


if __name__ == "__main__":
    main()
