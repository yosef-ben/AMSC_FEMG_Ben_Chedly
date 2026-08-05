#!/usr/bin/env python3

"""Check that every number shown in a report figure comes from a stored result.

Each check re-derives, from the CSV or VTK file the figure is drawn from, a
quantity that the figure states in text, and compares it with the value written
in the figure or in its README. The point is not to re-run the simulations but
to make sure that nothing displayed was typed by hand or produced by a
different run than the one stored.

Run from the project root; a non-zero exit status means a figure and its data
disagree.
"""

import csv
import glob
import json
import math
import re
import sys
from pathlib import Path

BENCH = Path("benchmarks")
TOLERANCE = 5e-3


def read_csv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def close(found, expected, tolerance=TOLERANCE):
    if expected == 0.0:
        return abs(found) <= tolerance
    return abs(found - expected) <= tolerance * max(1.0, abs(expected))


class Report:
    def __init__(self):
        self.rows = []
        self.failed = 0

    def check(self, figure, quantity, found, expected, unit=""):
        ok = close(found, expected)
        self.failed += not ok
        self.rows.append((figure, quantity, found, expected, unit, ok))

    def note(self, figure, quantity, message):
        self.rows.append((figure, quantity, message, "", "", True))

    def show(self):
        width = max(len(row[0]) for row in self.rows)
        for figure, quantity, found, expected, unit, ok in self.rows:
            mark = "ok  " if ok else "FAIL"
            if expected == "":
                print(f"{mark} {figure:<{width}}  {quantity}: {found}")
            else:
                print(f"{mark} {figure:<{width}}  {quantity}: "
                      f"{found:.6g} vs {expected:.6g} {unit}".rstrip())
        print(f"\n{len(self.rows)} checks, {self.failed} failing")


def crossing(times, values, level):
    for k in range(1, len(values)):
        if values[k - 1] < level <= values[k]:
            span = values[k] - values[k - 1]
            return times[k - 1] + (level - values[k - 1]) / span * (
                times[k] - times[k - 1])
    return float("nan")


def check_18(report):
    """Front speeds annotated on the 1D sensitivity figure."""
    name = "18 front_speeds"
    rows = read_csv(BENCH / "18_fisher_kolmogorov_1d_sensitivity/results"
                           "/front_speeds.csv")
    for row in rows:
        diffusion = float(row["diffusion"])
        alpha = float(row["alpha"])
        theory = 2.0 * math.sqrt(diffusion * alpha)
        report.check(name, f"2*sqrt(d*alpha) at d={diffusion:g}",
                     float(row["theoretical_speed"]), theory)
        stated = float(row["relative_error"])
        measured = float(row["measured_speed"])
        report.check(name, f"relative error at d={diffusion:g}",
                     (measured - theory) / theory, stated)

    # The time-step study quotes three front positions in its README.
    name = "18 time_step_study"
    rows = read_csv(BENCH / "18_fisher_kolmogorov_1d_sensitivity/results"
                           "/time_step_study.csv")
    positions = {float(row["dt"]): float(row["front_position"])
                 for row in rows}
    for step, expected in ((0.025, 0.6889), (0.05, 0.7268), (0.1, 0.8145)):
        report.check(name, f"front position at dt={step:g}", positions[step],
                     expected)


def check_19(report):
    """Lobe separation and refinement rates of the Fornari comparison."""
    name = "19 biomarker_comparison"
    base = BENCH / "19_fisher_kolmogorov_fornari83/results"
    for label, filename, expected in (
            ("nodal", "nodal_biomarkers.csv", 2.9469e-07),
            ("FEM", "fem_biomarkers.csv", 5.7e-03)):
        rows = read_csv(base / filename)
        times = [float(row["time"]) for row in rows]
        crossings = [crossing(times, [float(row[lobe]) for row in rows], 50.0)
                     for lobe in ("temporal", "frontal", "parietal",
                                  "occipital")]
        spread = max(crossings) - min(crossings)
        report.check(name, f"{label} lobe spread", spread, expected,
                     "years")

    name = "19 refinement"
    rows = read_csv(base / "space_refinement.csv")
    for row, expected in zip(rows, (0.4389, 0.1178, 0.0229, 0.0)):
        report.check(name, f"max difference, {row['cells_per_edge']} cells",
                     float(row["max_biomarker_difference"]), expected, "%")
    rows = read_csv(base / "time_refinement.csv")
    rates = [float(row["rate"]) for row in rows if row["rate"]]
    for rate, expected in zip(rates, (1.0552, 1.0263, 1.0126)):
        report.check(name, "observed temporal rate", rate, expected)


def check_20(report):
    """Crossing times of the conversion-rate study."""
    name = "20 alpha_sensitivity"
    rows = read_csv(BENCH / "20_fisher_kolmogorov_alpha_sensitivity/results"
                           "/alpha_sensitivity.csv")
    expected = {("Nodal reference", 0.5): 11.1054737623,
                ("Nodal reference", 0.1): 59.3302643285,
                ("Metric-graph FEM", 0.5): 12.6760711624,
                ("Metric-graph FEM", 0.1): 67.8056471471}
    for row in rows:
        key = (row["method"], float(row["alpha"]))
        if key in expected:
            report.check(name, f"t50 {key[0]} alpha={key[1]:g}",
                         float(row["time_50_percent"]), expected[key],
                         "years")
        if float(row["alpha"]) == 0.0:
            drift = abs(float(row["global_drift_alpha_zero"]))
            report.check(name, f"conservation drift {row['method']}",
                         drift, 0.0, "")


def check_21(report):
    """Regional means and their ranking against the published table."""
    name = "21 regional_averages"
    rows = read_csv(BENCH / "21_fisher_kolmogorov_corti83/results"
                           "/regional_averages.csv")
    final = rows[-1]
    report.check(name, "final time", float(final["time"]), 20.0, "years")
    expected = {"frontal": 0.202084, "temporal": 0.166959,
                "parietal": 0.096729, "insular": 0.129195,
                "limbic": 0.210109, "occipital": 0.073457,
                "subcortical": 0.127043}
    for region, value in expected.items():
        report.check(name, f"{region} at T", float(final[region]), value)

    corti = {"frontal": 0.9289, "limbic": 0.8905, "temporal": 0.8699,
             "insular": 0.8558, "subcortical": 0.8413, "parietal": 0.7738,
             "occipital": 0.7336}
    ours = {region: float(final[region]) for region in corti}
    order_corti = sorted(corti, key=lambda r: -corti[r])
    order_ours = sorted(ours, key=lambda r: -ours[r])
    rank_corti = {r: i for i, r in enumerate(order_corti)}
    rank_ours = {r: i for i, r in enumerate(order_ours)}
    concordant = sum(1 for a in corti for b in corti
                     if rank_corti[a] < rank_corti[b]
                     and rank_ours[a] < rank_ours[b])
    report.check(name, "concordant pairs vs Corti table 3", concordant, 20)
    squared = sum((rank_corti[r] - rank_ours[r]) ** 2 for r in corti)
    spearman = 1 - 6 * squared / (len(corti) * (len(corti) ** 2 - 1))
    report.check(name, "rank correlation", spearman, 0.9643)


def check_22(report):
    """Degrees of freedom, nonzeros and per-step cost of the baseline."""
    name = "22 sequential_performance"
    rows = read_csv(BENCH / "22_fisher_kolmogorov_sequential_performance"
                           "/results/sequential_performance.csv")
    expected = {1: (83, 2343), 2: (1213, 5733), 4: (3473, 12513),
                8: (7993, 26073)}
    for row in rows:
        cells = int(row["cells_per_edge"])
        dofs, nonzeros = expected[cells]
        report.check(name, f"DoFs, {cells} cells", int(row["n_dofs"]), dofs)
        report.check(name, f"nonzeros, {cells} cells", int(row["matrix_nnz"]),
                     nonzeros)
        report.check(name, f"seconds per step, {cells} cells",
                     float(row["seconds_per_step"]),
                     float(row["solve_seconds"]) / 100.0)


def check_23(report):
    """Damkohler numbers and the boundedness limit of the FEM."""
    name = "23 diffusion_scaling"
    base = BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
    with open(base / "diffusion_scaling_summary.json") as stream:
        summary = json.load(stream)
    report.check(name, "Fiedler value", summary["fiedler_value"], 0.795445)
    report.check(name, "benchmark 21 scaling", summary["benchmark_21_scaling"],
                 1.0 / summary["maximum_adjacency"])
    rows = read_csv(base / "diffusion_scaling.csv")
    for row in rows:
        scaling = float(row["diffusion_scaling"])
        report.check(name, f"Da at rho={scaling:g}", float(row["damkohler"]),
                     summary["alpha"] / (scaling * summary["fiedler_value"]))
    bounded = [float(row["diffusion_scaling"]) for row in rows
               if row["fem_corti_crank_nicolson"] == "bounded"]
    report.check(name, "smallest bounded scaling", min(bounded), 0.05)


def check_25(report):
    """Extremes and spans of the seeding study."""
    name = "25 seeding_vulnerability"
    base = BENCH / "25_connectome_seeding_vulnerability/results"
    for scaling, span, fastest, slowest in (
            (1, 0.018, "Right-Caudate", "ctx-rh-frontalpole"),
            (0.02, 4.020, "Right-Caudate", "ctx-rh-temporalpole")):
        rows = read_csv(base / f"seeding_vulnerability_rho_{scaling:g}.csv")
        values = [float(row["infection_time_years"]) for row in rows]
        report.check(name, f"span at rho={scaling:g}",
                     max(values) - min(values), span, "years")
        order = sorted(rows, key=lambda row: float(
            row["infection_time_years"]))
        report.note(name, f"fastest at rho={scaling:g}",
                    f"{order[0]['name']} (expected {fastest})")
        report.note(name, f"slowest at rho={scaling:g}",
                    f"{order[-1]['name']} (expected {slowest})")
        if order[0]["name"] != fastest or order[-1]["name"] != slowest:
            report.failed += 1
        entorhinal = [index + 1 for index, row in enumerate(order)
                      if "entorhinal" in row["name"]]
        report.note(name, f"entorhinal ranks at rho={scaling:g}",
                    f"{entorhinal} of {len(order)}")


def check_26(report):
    """Field range of the anatomical progression."""
    name = "26 anatomical_progression"
    try:
        import vtk
    except ImportError:
        report.note(name, "field range", "skipped, vtk not installed")
        return
    paths = sorted(glob.glob(
        "output/fisher_kolmogorov/corti83_refined/solution_*.vtp"))
    if not paths:
        report.note(name, "field range", "skipped, run the benchmark first")
        return
    low, high, points = math.inf, -math.inf, 0
    for path in paths:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(path)
        reader.Update()
        array = reader.GetOutput().GetPointData().GetArray("c")
        low = min(low, array.GetRange()[0])
        high = max(high, array.GetRange()[1])
        points = reader.GetOutput().GetNumberOfPoints()
    report.check(name, "sampling points", points, 10170)
    report.check(name, "field minimum", low, 0.01)
    report.check(name, "field maximum", high, 0.9987)
    report.check(name, "stays within [0,1]", 1.0 if 0.0 <= low and high <= 1.0
                 else 0.0, 1.0)

    # The three stages are selected by a rule, not by eye: the stored times at
    # which the network mean reaches 0.2, 0.4 and 0.6.
    import numpy as np
    means = {}
    for path in paths:
        step = int(re.search(r"solution_(\d+)\.vtp$", path).group(1))
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(path)
        reader.Update()
        array = reader.GetOutput().GetPointData().GetArray("c")
        means[step] = float(np.mean([array.GetValue(i)
                                     for i in range(array.GetNumberOfTuples())]))
    for target, expected in ((0.2, 23.0), (0.4, 31.0), (0.6, 38.0)):
        step = min(means, key=lambda s: abs(means[s] - target))
        report.check(name, f"stage at mean c = {target:g}", step * 0.2,
                     expected, "years")


def main():
    report = Report()
    for check in (check_18, check_19, check_20, check_21, check_22, check_23,
                  check_25, check_26):
        try:
            check(report)
        except FileNotFoundError as error:
            report.note(check.__name__, "input", f"missing: {error.filename}")
            report.failed += 1
    report.show()
    return 1 if report.failed else 0


if __name__ == "__main__":
    sys.exit(main())
