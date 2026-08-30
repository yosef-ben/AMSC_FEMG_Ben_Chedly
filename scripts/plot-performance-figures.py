#!/usr/bin/env python3

"""The three figures of the performance chapter, house style.

Ordering: factor fill and per-step cost of the candidate orderings against
the unknown count, from the sequential ordering study. Condensation: the
condensed step against the full-system reference and the phase split of the
condensed step, from the static-condensation study. Scaling: OpenMP and MPI
speedup on the 128-cell connectome and the equal-worker hybrid comparison.
Every number is read from the stored benchmark CSV records; the derived
values quoted in the chapter text are printed for cross-checking.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style

ORDERING = Path("benchmarks/28_sequential_ordering_study/results"
                "/ordering_study.csv")
CONDENSATION = Path("benchmarks/29_static_condensation/results"
                    "/condensation_study.csv")
OPENMP = Path("benchmarks/30_openmp_condensation/results/omp_scaling.csv")
MPI = Path("benchmarks/31_mpi_condensation/results/mpi_scaling.csv")
HYBRID = Path("benchmarks/32_hybrid_condensation/results"
              "/hybrid_comparison.csv")

VARIANTS = (
    ("lu_colamd", "LU, COLAMD", "#D62728", "o"),
    ("lu_natural", "LU, natural", "#FF7F0E", "s"),
    ("lu_rcm", "LU, RCM", "#9467BD", "^"),
    ("ldlt_natural", "LDLT, natural", "#1F77B4", "D"),
)


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def save(figure, output_dir, stem):
    for suffix in (".pdf", ".png"):
        figure.savefig(output_dir / f"{stem}{suffix}", dpi=240,
                       bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"Saved {output_dir / stem}.pdf")


def ordering_figure(output_dir):
    rows = read_rows(ORDERING)
    series = {}
    for row in rows:
        series.setdefault(row["variant"], []).append(row)

    figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    shifts_fill = {"lu_colamd": 0.62, "lu_natural": 1.55, "lu_rcm": 1.0,
                   "ldlt_natural": 0.72}
    shifts_time = {"lu_colamd": 1.35, "lu_natural": 0.72, "lu_rcm": 1.0,
                   "ldlt_natural": 0.85}
    for key, label, colour, marker in VARIANTS:
        data = series[key]
        dofs = [float(row["n_dofs"]) for row in data]
        fill = [float(row["factor_nnz"]) for row in data]
        step = [float(row["step_median_s"]) for row in data]
        axes[0].loglog(dofs, fill, marker=marker, linestyle="-",
                       color=colour, linewidth=1.8, markersize=5.0,
                       markeredgecolor="white", markeredgewidth=0.7)
        figure_style.label_series(axes[0], dofs[-1] * 1.3,
                                  fill[-1] * shifts_fill[key], label,
                                  colour, fontsize=8.5)
        axes[1].loglog(dofs, step, marker=marker, linestyle="-",
                       color=colour, linewidth=1.8, markersize=5.0,
                       markeredgecolor="white", markeredgewidth=0.7)
        figure_style.label_series(axes[1], dofs[-1] * 1.3,
                                  step[-1] * shifts_time[key], label,
                                  colour, fontsize=8.5)

    colamd = {int(row["n_dofs"]): float(row["factor_nnz"])
              for row in series["lu_colamd"]}
    # The top-left corner is the only region the rising curves never
    # reach, so the note lives there and the arrow drops onto the bump.
    axes[0].annotate("2 elements/connection:\nmore fill than the\n"
                     "next two meshes",
                     xy=(1213, colamd[1213]), xytext=(70, 2.2e7),
                     va="top", fontsize=8, fontweight="bold", color="0.35",
                     arrowprops={"arrowstyle": "->", "color": "0.35",
                                 "linewidth": 1.1})
    axes[0].set_xlim(45, 2.2e6)
    axes[0].set_ylim(8e2, 6e7)
    axes[0].set_ylabel("factor nonzeros", labelpad=2)
    axes[1].set_xlim(45, 2.2e6)
    axes[1].set_ylim(8e-5, 30)
    axes[1].set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    axes[1].set_yticklabels(["0.0001", "0.001", "0.01", "0.1", "1"])
    axes[1].set_ylabel("seconds per time step", labelpad=2)
    for axis in axes:
        axis.minorticks_off()
        figure_style.xname(axis, "degrees of freedom", y=-0.12)
    for letter, axis in zip("ab", axes):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")
    figure.tight_layout(w_pad=3.2)
    save(figure, output_dir, "performance_ordering")

    print("  chapter table (cells, variant, factor nnz, step ms):")
    for cells in ("2", "8", "64"):
        for key, label, _, _ in VARIANTS:
            row = next(r for r in series[key]
                       if r["cells_per_edge"] == cells)
            print(f"    {cells:>2} {label:<14} {row['n_dofs']:>6} "
                  f"{int(row['factor_nnz']):>9} "
                  f"{float(row['step_median_s'])*1e3:8.2f}")
    at64 = {key: next(r for r in series[key] if r["cells_per_edge"] == "64")
            for key in series}
    print(f"  colamd fill 1213 vs 3473: {colamd[1213]:.0f} vs "
          f"{colamd[3473]:.0f}")
    print(f"  rcm/natural fill at 64: "
          f"{float(at64['lu_rcm']['factor_nnz']) / float(at64['lu_natural']['factor_nnz']):.1f}")
    print(f"  rcm/natural profile at 64: "
          f"{float(at64['lu_rcm']['profile']) / float(at64['lu_natural']['profile']):.1f}")
    print(f"  colamd/ldlt_natural step at 64: "
          f"{float(at64['lu_colamd']['step_median_s']) / float(at64['ldlt_natural']['step_median_s']):.2f}")
    ldlt = at64["ldlt_natural"]
    print(f"  ldlt_natural at 64: rebuild "
          f"{float(ldlt['rebuild_median_s'])*1e3:.1f} factor "
          f"{float(ldlt['factor_median_s'])*1e3:.1f} solve "
          f"{float(ldlt['solve_median_s'])*1e3:.2f} ms")


def condensation_figure(output_dir):
    rows = read_rows(CONDENSATION)
    dofs = [float(row["n_dofs"]) for row in rows]
    condensed = [float(row["cond_step_median_s"]) for row in rows]
    full = [float(row["full_step_median_s"]) for row in rows]

    figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    axis = axes[0]
    axis.loglog(dofs, full, marker="o", linestyle="-", color="#1F77B4",
                linewidth=1.8, markersize=5.5, markeredgecolor="white",
                markeredgewidth=0.7)
    axis.loglog(dofs, condensed, marker="s", linestyle="-", color="#D62728",
                linewidth=1.8, markersize=5.5, markeredgecolor="white",
                markeredgewidth=0.7)
    figure_style.label_series(axis, dofs[-1] * 1.3, full[-1],
                              "full system,\nLDLT natural", "#1F77B4",
                              fontsize=8.5)
    figure_style.label_series(axis, dofs[-1] * 1.3, condensed[-1],
                              "condensed", "#D62728", fontsize=8.5)
    axis.annotate("", xy=(dofs[-1], full[-1] * 0.82),
                  xytext=(dofs[-1], condensed[-1] * 1.22),
                  arrowprops={"arrowstyle": "<->", "color": "0.35",
                              "linewidth": 1.1})
    axis.text(dofs[-1] * 0.82, (full[-1] * condensed[-1]) ** 0.5,
              f"{full[-1] / condensed[-1]:.1f}x", fontsize=8.5,
              fontweight="bold", color="0.35", ha="right", va="center")
    axis.set_xlim(4.5e3, 1.5e6)
    axis.set_ylim(5e-4, 0.3)
    axis.set_yticks([1e-3, 1e-2, 1e-1])
    axis.set_yticklabels(["0.001", "0.01", "0.1"])
    axis.set_ylabel("seconds per time step", labelpad=2)

    axis = axes[1]
    phases = (
        ("cond_local_median_s", "per-connection\nelimination", "#D62728",
         "o", 1.0),
        ("cond_back_median_s", "back substitution", "#FF7F0E", "s", 1.0),
        ("cond_interface_median_s", "interface solve", "#1F77B4", "D", 1.0),
    )
    for name, label, colour, marker, shift in phases:
        values = [float(row[name]) for row in rows]
        axis.loglog(dofs, values, marker=marker, linestyle="-",
                    color=colour, linewidth=1.8, markersize=5.0,
                    markeredgecolor="white", markeredgewidth=0.7)
        figure_style.label_series(axis, dofs[-1] * 1.3, values[-1] * shift,
                                  label, colour, fontsize=8.5)
    axis.set_xlim(4.5e3, 1.5e6)
    axis.set_ylim(8e-6, 5e-2)
    axis.set_yticks([1e-5, 1e-4, 1e-3, 1e-2])
    axis.set_yticklabels(["0.00001", "0.0001", "0.001", "0.01"])
    axis.set_ylabel("seconds per phase", labelpad=2)

    for axis in axes:
        axis.minorticks_off()
        figure_style.xname(axis, "degrees of freedom", y=-0.12)
    for letter, axis in zip("ab", axes):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")
    figure.tight_layout(w_pad=3.2)
    save(figure, output_dir, "performance_condensation")

    print("  chapter table (cells, n, full ms, condensed ms, ratio):")
    for row in rows:
        print(f"    {row['cells_per_edge']:>3} {row['n_dofs']:>6} "
              f"{float(row['full_step_median_s'])*1e3:6.2f} "
              f"{float(row['cond_step_median_s'])*1e3:6.2f} "
              f"{float(row['full_step_median_s']) / float(row['cond_step_median_s']):.1f}")
    last = rows[-1]
    print(f"  phases at 128: local "
          f"{float(last['cond_local_median_s'])*1e3:.2f} interface "
          f"{float(last['cond_interface_median_s'])*1e3:.3f} back "
          f"{float(last['cond_back_median_s'])*1e3:.2f} ms, interface "
          f"share {float(last['cond_interface_median_s']) / float(last['cond_step_median_s']) * 100:.1f}%")


def scaling_figure(output_dir):
    omp = read_rows(OPENMP)
    mpi = read_rows(MPI)
    hybrid = read_rows(HYBRID)

    figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    axis = axes[0]
    workers = [1, 2, 4, 8]
    axis.plot(workers, workers, linestyle="--", color="0.6", linewidth=1.3)
    figure_style.label_series(axis, 6.7, 7.35, "ideal", "0.5", fontsize=8.5)
    axis.axvline(4.0, color="0.85", linewidth=1.1, zorder=0)
    axis.text(4.12, 0.22, "physical cores", fontsize=8, fontweight="bold",
              color="0.45", ha="left", va="bottom")
    for rows, count_key, step_key, label, colour, marker, shift in (
            (omp, "threads", "omp_step_median_s", "OpenMP", "#D62728",
             "o", 0.35),
            (mpi, "ranks", "step_median_s", "MPI", "#1F77B4", "s", -0.5)):
        points = [(int(row[count_key]), float(row["speedup"]))
                  for row in rows if int(row[count_key]) > 1]
        axis.plot([p[0] for p in points], [p[1] for p in points],
                  marker=marker, linestyle="-", color=colour,
                  linewidth=1.8, markersize=5.5, markeredgecolor="white",
                  markeredgewidth=0.7)
        figure_style.label_series(axis, 8.25, points[-1][1] + shift, label,
                                  colour, fontsize=9)
    axis.set_xlim(0.6, 9.4)
    axis.set_ylim(0.0, 8.6)
    axis.set_xticks(workers)
    axis.set_yticks([0, 2, 4, 6, 8])
    axis.set_ylabel("speedup", labelpad=2)
    figure_style.xname(axis, "workers", y=-0.12)

    axis = axes[1]
    labels = [f"{row['ranks']}x{row['threads']}" for row in hybrid]
    times = [float(row["step_median_s"]) * 1e3 for row in hybrid]
    bars = axis.bar(range(len(hybrid)), times, width=0.62,
                    color=["#D62728", "#B04A58", "#7A5D87", "#1F77B4"])
    for bar, time, row in zip(bars, times, hybrid):
        axis.text(bar.get_x() + bar.get_width() / 2.0, time + 0.06,
                  f"{time:.2f} ms\nS = {float(row['speedup']):.2f}",
                  fontsize=8.5, fontweight="bold", color="0.25",
                  ha="center", va="bottom", linespacing=1.25)
    axis.set_xticks(range(len(hybrid)))
    axis.set_xticklabels(labels)
    axis.set_ylim(0.0, 3.1)
    axis.set_ylabel("milliseconds per time step", labelpad=2)
    figure_style.xname(axis, "ranks x threads", y=-0.12)

    for letter, axis in zip("ab", axes):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")
    figure.tight_layout(w_pad=3.2)
    save(figure, output_dir, "performance_scaling")

    for name, rows, count_key, step_key in (
            ("OpenMP", omp, "threads", "omp_step_median_s"),
            ("MPI", mpi, "ranks", "step_median_s")):
        t1 = float(rows[0]["t1_pooled_s"])
        print(f"  {name}: sequential reference {t1*1e3:.2f} ms")
        for row in rows:
            print(f"    {count_key} {row[count_key]}: "
                  f"{float(row[step_key])*1e3:5.2f} ms  "
                  f"S {float(row['speedup']):.2f}  "
                  f"E {float(row['efficiency']):.2f}")
    t1 = float(hybrid[0]["t1_pooled_s"])
    print(f"  hybrid: sequential reference {t1*1e3:.2f} ms")
    for row in hybrid:
        print(f"    {row['ranks']}x{row['threads']}: "
              f"{float(row['step_median_s'])*1e3:.2f} ms  "
              f"S {float(row['speedup']):.2f}  allreduce "
              f"{float(row['allreduce_median_s'])*1e3:.2f} ms")
    omp8 = next(row for row in omp if row["threads"] == "8")
    print(f"  OpenMP 8 phases: local {float(omp8['omp_local_median_s'])*1e3:.2f} "
          f"reduce {float(omp8['omp_reduce_median_s'])*1e3:.3f} "
          f"iface {float(omp8['omp_interface_median_s'])*1e3:.3f} "
          f"back {float(omp8['omp_back_median_s'])*1e3:.3f} ms")
    omp1 = next(row for row in omp if row["threads"] == "1")
    print(f"  OpenMP 1-thread overhead: "
          f"{(float(omp1['omp_step_median_s']) / float(omp1['t1_pooled_s']) - 1.0) * 100:.1f}%")
    mpi1 = next(row for row in mpi if row["ranks"] == "1")
    print(f"  MPI 1-rank overhead: "
          f"{(float(mpi1['step_median_s']) / float(mpi1['t1_pooled_s']) - 1.0) * 100:.1f}%")
    mpi4 = next(row for row in mpi if row["ranks"] == "4")
    print(f"  MPI allreduce at 2/4/8 ranks: "
          + "/".join(f"{float(row['allreduce_median_s'])*1e3:.2f}"
                     for row in mpi if row["ranks"] != "1") + " ms")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("report"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_style.apply()
    ordering_figure(args.output_dir)
    condensation_figure(args.output_dir)
    scaling_figure(args.output_dir)


if __name__ == "__main__":
    main()
