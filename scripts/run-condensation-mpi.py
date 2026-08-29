#!/usr/bin/env python3

"""Strong scaling of the MPI condensation on the 128-cell connectome.

Ranks own contiguous edge blocks balanced by interior unknowns; the 83
vertices are the replicated interface, summed by two MPI_Allreduce calls
per step. Up to four ranks each process is bound to its own physical core;
eight ranks use the hardware threads, matching the OpenMP study. The
sequential reference T1 is pooled from repeated runs of the sequential
condensed engine (test_condensation) in the same session, and every field
is the median across three repeated processes.
"""

import argparse
import csv
import os
import statistics
import subprocess
from pathlib import Path

RAW_FIELDS = [
    "n_dofs", "ranks", "steps",
    "local_median_s", "allreduce_median_s", "interface_median_s",
    "back_median_s", "step_median_s", "loop_total_s",
]
RANKS = (1, 2, 4, 8)
VALIDATIONS = ((8, 2), (8, 4), (8, 8), (128, 8))
SCALING_CELLS = 128
REPEATS = 3


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build-release/test_condensation_mpi"))
    parser.add_argument("--sequential", type=Path,
                        default=Path("build-release/test_condensation"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=200)
    return parser.parse_args()


def environment():
    env = os.environ.copy()
    env.update({
        "OMPI_ALLOW_RUN_AS_ROOT": "1",
        "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    return env


def run_mpi(executable, cells, mode, ranks, max_steps):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    binding = (["--map-by", "core", "--bind-to", "core"] if ranks <= 4
               else ["--use-hwthread-cpus", "--map-by", "hwthread",
                     "--bind-to", "hwthread"])
    completed = subprocess.run(
        ["mpiexec", *binding, "-n", str(ranks), str(executable),
         str(graph), mode, str(max_steps)],
        check=True, text=True, capture_output=True, env=environment(),
    )
    return completed.stdout


def run_sequential(executable, cells, max_steps):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    completed = subprocess.run(
        [str(executable), str(graph), "benchmark", str(max_steps)],
        check=True, text=True, capture_output=True, env=environment(),
    )
    line = next(line for line in completed.stdout.splitlines()
                if line.startswith("CONDENSE,"))
    values = line.split(",")[1:]
    # cond_step_median is the sum of the three condensed phase medians.
    return float(values[5])


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.output_dir / "mpi_validation.csv", "w",
              newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", "n_dofs", "ranks", "steps",
                         "max_diff_mpi_vs_seq"])
        for cells, ranks in VALIDATIONS:
            stdout = run_mpi(args.executable, cells, "--validate", ranks,
                             100)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("VALIDATION,"))
            writer.writerow([cells, *line.split(",")[1:]])
            print(f"validate cells={cells} ranks={ranks}: "
                  f"max diff {line.split(',')[-1]}")

    t1_samples = [run_sequential(args.sequential, SCALING_CELLS,
                                 args.max_steps) for _ in range(REPEATS)]
    t1 = statistics.median(t1_samples)
    print(f"pooled sequential T1: {t1*1e3:.2f} ms per step")

    with open(args.output_dir / "mpi_scaling.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", *RAW_FIELDS,
                         "t1_pooled_s", "speedup", "efficiency"])
        for ranks in RANKS:
            repeats = []
            for _ in range(REPEATS):
                stdout = run_mpi(args.executable, SCALING_CELLS,
                                 "benchmark", ranks, args.max_steps)
                line = next(line for line in stdout.splitlines()
                            if line.startswith("MPI,"))
                repeats.append(dict(zip(RAW_FIELDS,
                                        line.split(",")[1:])))
            record = {}
            for name in RAW_FIELDS:
                if name in ("n_dofs", "ranks", "steps"):
                    record[name] = repeats[0][name]
                else:
                    record[name] = statistics.median(
                        float(r[name]) for r in repeats)
            speedup = t1 / float(record["step_median_s"])
            efficiency = speedup / ranks
            writer.writerow([SCALING_CELLS,
                             *[record[name] for name in RAW_FIELDS],
                             f"{t1:.6f}",
                             f"{speedup:.4f}", f"{efficiency:.4f}"])
            print(f"ranks={ranks}: "
                  f"step {float(record['step_median_s'])*1e3:6.2f} ms "
                  f"(local {float(record['local_median_s'])*1e3:6.2f}"
                  f", allreduce {float(record['allreduce_median_s'])*1e3:5.2f}"
                  f", iface {float(record['interface_median_s'])*1e3:5.2f}"
                  f", back {float(record['back_median_s'])*1e3:5.2f})  "
                  f"S={speedup:.2f} E={efficiency:.2f}")


if __name__ == "__main__":
    main()
