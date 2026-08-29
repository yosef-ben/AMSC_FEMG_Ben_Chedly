#!/usr/bin/env python3

"""Hybrid MPI+OpenMP comparison at equal total worker count.

Four configurations of the same eight hardware threads of the 4-core
machine, measured with the one hybrid binary in one session: pure OpenMP
(1x8), two hybrids (2x4, 4x2) and pure MPI (8x1). Each rank is packed on
its own block of hardware threads (--map-by slot:PE=threads) and OpenMP
threads bind inside it. The sequential reference T1 is pooled from
repeated runs of the sequential condensed engine.
"""

import argparse
import csv
import os
import statistics
import subprocess
from pathlib import Path

RAW_FIELDS = [
    "n_dofs", "ranks", "threads", "steps",
    "local_median_s", "allreduce_median_s", "interface_median_s",
    "back_median_s", "step_median_s", "loop_total_s",
]
CONFIGS = ((1, 8), (2, 4), (4, 2), (8, 1))
SCALING_CELLS = 128
REPEATS = 3


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build-release/test_condensation_hybrid"))
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
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    return env


def run_hybrid(executable, cells, mode, ranks, threads, max_steps):
    graph = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt")
    completed = subprocess.run(
        ["mpiexec", "--use-hwthread-cpus",
         "--map-by", f"slot:PE={threads}", "--bind-to", "hwthread",
         "-n", str(ranks),
         "-x", f"OMP_NUM_THREADS={threads}",
         "-x", "OMP_PROC_BIND=close", "-x", "OMP_PLACES=threads",
         str(executable), str(graph), mode, str(max_steps), str(threads)],
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
    return float(line.split(",")[1:][5])


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.output_dir / "hybrid_validation.csv", "w",
              newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", "n_dofs", "ranks", "threads",
                         "steps", "max_diff_hybrid_vs_seq"])
        for cells, (ranks, threads) in (
                (8, (2, 4)), (8, (4, 2)), (128, (2, 4))):
            stdout = run_hybrid(args.executable, cells, "--validate",
                                ranks, threads, 100)
            line = next(line for line in stdout.splitlines()
                        if line.startswith("VALIDATION,"))
            writer.writerow([cells, *line.split(",")[1:]])
            print(f"validate cells={cells} {ranks}x{threads}: "
                  f"max diff {line.split(',')[-1]}")

    t1 = statistics.median(
        [run_sequential(args.sequential, SCALING_CELLS, args.max_steps)
         for _ in range(REPEATS)])
    print(f"pooled sequential T1: {t1*1e3:.2f} ms per step")

    with open(args.output_dir / "hybrid_comparison.csv", "w",
              newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cells_per_edge", *RAW_FIELDS,
                         "t1_pooled_s", "speedup", "efficiency"])
        for ranks, threads in CONFIGS:
            repeats = []
            for _ in range(REPEATS):
                stdout = run_hybrid(args.executable, SCALING_CELLS,
                                    "benchmark", ranks, threads,
                                    args.max_steps)
                line = next(line for line in stdout.splitlines()
                            if line.startswith("HYBRID,"))
                repeats.append(dict(zip(RAW_FIELDS,
                                        line.split(",")[1:])))
            record = {}
            for name in RAW_FIELDS:
                if name in ("n_dofs", "ranks", "threads", "steps"):
                    record[name] = repeats[0][name]
                else:
                    record[name] = statistics.median(
                        float(r[name]) for r in repeats)
            speedup = t1 / float(record["step_median_s"])
            efficiency = speedup / 8.0
            writer.writerow([SCALING_CELLS,
                             *[record[name] for name in RAW_FIELDS],
                             f"{t1:.6f}",
                             f"{speedup:.4f}", f"{efficiency:.4f}"])
            print(f"{ranks}x{threads}: "
                  f"step {float(record['step_median_s'])*1e3:6.2f} ms "
                  f"(local {float(record['local_median_s'])*1e3:6.2f}"
                  f", allreduce {float(record['allreduce_median_s'])*1e3:5.2f}"
                  f", iface {float(record['interface_median_s'])*1e3:5.2f}"
                  f", back {float(record['back_median_s'])*1e3:5.2f})  "
                  f"S={speedup:.2f} E={efficiency:.2f}")


if __name__ == "__main__":
    main()
