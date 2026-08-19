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

    def check_contains(self, figure, quantity, text, needle):
        ok = needle in text
        self.failed += not ok
        self.rows.append((figure, quantity,
                          f"'{needle}' {'found' if ok else 'MISSING'}",
                          "", "", ok))

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


def check_19_topology(report):
    """Ranges and structure the connectome topology figure displays."""
    name = "19 connectome_topology"
    import numpy as np
    edges = read_csv(Path("data/connectome/fornari83/edges.csv"))
    size = 83
    adjacency = np.zeros((size, size))
    for row in edges:
        i, j = int(row["source"]), int(row["target"])
        adjacency[i, j] = adjacency[j, i] = float(row["connectivity_weight"])
    degree = (adjacency > 0).sum(axis=1)
    weighted = adjacency.sum(axis=1)
    nonzero = adjacency[adjacency > 0]

    # The three ranges printed at the ends of the two colour bars, which are
    # also the values Fornari et al. publish: 6-48, 2.1-127.6 and 0.01-35.32.
    report.check(name, "smallest degree", float(degree.min()), 6.0)
    report.check(name, "largest degree", float(degree.max()), 48.0)
    report.check(name, "smallest weighted degree", float(weighted.min()),
                 2.05052)
    report.check(name, "largest weighted degree", float(weighted.max()),
                 127.6435)
    report.check(name, "mean weighted degree", float(weighted.mean()),
                 42.7547)
    report.check(name, "smallest adjacency", float(nonzero.min()), 0.0084602)
    report.check(name, "largest adjacency", float(nonzero.max()), 35.3221)
    report.check(name, "mean adjacency", float(nonzero.mean()), 1.57019)

    report.check(name, "symmetric", 1.0 if np.allclose(adjacency,
                                                       adjacency.T) else 0.0,
                 1.0)
    report.check(name, "no self-loops",
                 1.0 if np.all(np.diag(adjacency) == 0) else 0.0, 1.0)
    report.check(name, "non-zero cells", float(len(nonzero)), 2260.0)

    # The regions attaining the extremes, which the reference names: the
    # frontal pole and the precentral gyrus for the weighted degree, the
    # superior parietal to precuneus pair for the largest adjacency and the
    # most fibres. The published pairing of the smallest and largest
    # adjacency is internally inconsistent with its own fibre counts; the
    # reconstruction resolves it, and these checks pin the resolution down.
    nodes = {int(row["node_id"]): row["name"]
             for row in read_csv(Path("data/connectome/fornari83/nodes.csv"))}
    identities = (
        ("smallest weighted degree region", nodes[int(weighted.argmin())],
         "ctx-rh-frontalpole"),
        ("largest weighted degree region", nodes[int(weighted.argmax())],
         "ctx-rh-precentral"),
    )
    masked = np.where(adjacency > 0, adjacency, np.inf)
    pair = np.unravel_index(np.argmax(adjacency), adjacency.shape)
    identities += (("largest adjacency pair",
                    " -- ".join(sorted(nodes[k] for k in pair)),
                    "ctx-rh-precuneus -- ctx-rh-superiorparietal"),)
    pair = np.unravel_index(np.argmin(masked), adjacency.shape)
    identities += (("smallest adjacency pair",
                    " -- ".join(sorted(nodes[k] for k in pair)),
                    "ctx-lh-isthmuscingulate -- ctx-lh-lateralorbitofrontal"),)
    fibre = np.zeros_like(adjacency)
    for row in edges:
        i, j = int(row["source"]), int(row["target"])
        fibre[i, j] = fibre[j, i] = float(row["fibre_number"])
    pair = np.unravel_index(np.argmax(fibre), fibre.shape)
    identities += (("largest fibre number pair",
                    " -- ".join(sorted(nodes[k] for k in pair)),
                    "ctx-rh-precuneus -- ctx-rh-superiorparietal"),)
    for label, found, expected in identities:
        report.note(name, label, f"{found} (expected {expected})")
        if found != expected:
            report.failed += 1
    report.check(name, "largest fibre number", float(fibre.max()), 595.5)

    # The colour scale is logarithmic; this is the measurement that justifies
    # it, and it is quoted in the caption.
    report.check(name, "decades spanned",
                 float(np.log10(nonzero.max() / nonzero.min())), 3.6207)
    report.check(name, "fraction below 5% of the maximum",
                 float((nonzero < 0.05 * nonzero.max()).mean()), 0.7735)

    # The eight warm clusters along the diagonal are the intra-lobe
    # connections; the report says their mean adjacency is between roughly
    # two and four times the mean of their hemisphere block.
    hemisphere = {int(row["node_id"]): row["hemisphere"]
                  for row in read_csv(
                      Path("data/connectome/fornari83/nodes.csv"))}
    group = {int(row["node_id"]): row["region"]
             for row in read_csv(BENCH / "21_fisher_kolmogorov_corti83"
                                         "/results/reaction_coefficients.csv")}
    ratios = []
    for side, low, high in (("right", 0, 41), ("left", 41, 82)):
        block = adjacency[low:high, low:high]
        block_mean = block[block > 0].mean()
        for lobe in ("frontal", "temporal", "parietal", "occipital"):
            members = [k for k in range(low, high)
                       if hemisphere[k] == side and group[k] == lobe]
            inside = adjacency[np.ix_(members, members)]
            ratios.append(float(inside[inside > 0].mean() / block_mean))
    report.check(name, "smallest intra-lobe to block mean ratio",
                 min(ratios), 1.912)
    report.check(name, "largest intra-lobe to block mean ratio",
                 max(ratios), 4.252)


def check_connectome_consistency(report):
    """The graph the solvers read is the graph the report describes.

    Recomputes, from nodes.csv and edges.csv, the statistics the connectome
    table quotes against the reference (degree, fibre number, fibre length,
    adjacency, weighted degree, with the regions at the extremes), the two
    structural fractions the text quotes, the Fiedler value the Damkohler
    section uses, and the correspondence between the weight file and the
    finite element meshes: same 1130 connections in the same order, unit
    metric length, and the stated number of cells. The longest retained fine
    connection is read from the summary written by the preparation script.
    """
    name = "connectome consistency"
    import numpy as np
    nodes = read_csv(Path("data/connectome/fornari83/nodes.csv"))
    edges = read_csv(Path("data/connectome/fornari83/edges.csv"))
    names = {int(row["node_id"]): row["name"] for row in nodes}
    hemisphere = {int(row["node_id"]): row["hemisphere"] for row in nodes}
    size = len(nodes)
    report.check(name, "vertices", float(size), 83.0)
    report.check(name, "edges", float(len(edges)), 1130.0)

    source = np.array([int(row["source"]) for row in edges])
    target = np.array([int(row["target"]) for row in edges])
    weight = np.array([float(row["connectivity_weight"]) for row in edges])
    number = np.array([float(row["fibre_number"]) for row in edges])
    length = np.array([float(row["fibre_length_mm"]) for row in edges])
    pairs = {(min(s, t), max(s, t)) for s, t in zip(source, target)}
    report.check(name, "simple graph (no loops, no repeated pairs)",
                 1.0 if len(pairs) == len(edges)
                 and np.all(source != target) else 0.0, 1.0)

    degree = np.zeros(size, int)
    weighted = np.zeros(size)
    for s, t, w in zip(source, target, weight):
        degree[s] += 1
        degree[t] += 1
        weighted[s] += w
        weighted[t] += w
    report.check(name, "degree range low", float(degree.min()), 6.0)
    report.check(name, "degree range high", float(degree.max()), 48.0)
    report.check(name, "fibre number mean", float(number.mean()), 40.1619)
    report.check(name, "fibre number range", float(number.max()), 595.5)
    report.check(name, "fibre length mean [mm]", float(length.mean()),
                 38.4009)
    report.check(name, "fibre length low [mm]", float(length.min()), 11.2867)
    report.check(name, "fibre length high, region pairs [mm]",
                 float(length.max()), 121.0235)
    for label, found, expected in (
            ("degree minimum region", names[int(degree.argmin())],
             "ctx-rh-frontalpole"),
            ("degree maximum region", names[int(degree.argmax())],
             "Right-Caudate"),
            ("longest mean-length pair",
             " -- ".join(sorted((names[int(source[length.argmax()])],
                                 names[int(target[length.argmax()])]))),
             "ctx-lh-lateralorbitofrontal -- ctx-lh-precuneus")):
        report.note(name, label, f"{found} (expected {expected})")
        if found != expected:
            report.failed += 1

    summary = json.loads(
        Path("data/connectome/fornari83/summary.json").read_text())
    report.check(name, "longest retained fine connection [mm]",
                 float(summary["fine_graph"]["fibre_length_mm"]["maximum"]),
                 136.8333)

    intra = sum(w for s, t, w in zip(source, target, weight)
                if hemisphere[s] == hemisphere[t]) / weight.sum()
    report.check(name, "weight fraction within a hemisphere", float(intra),
                 0.984)
    lobe_of = {int(row["node_id"]): row["region"]
               for row in read_csv(BENCH / "21_fisher_kolmogorov_corti83"
                                           "/results/reaction_coefficients.csv")}
    lobes = ("frontal", "temporal", "parietal", "occipital")
    inside = [w for s, t, w in zip(source, target, weight)
              if lobe_of[s] in lobes and lobe_of[s] == lobe_of[t]
              and hemisphere[s] == hemisphere[t]]
    outside = [w for s, t, w in zip(source, target, weight)
               if not (lobe_of[s] in lobes and lobe_of[s] == lobe_of[t]
                       and hemisphere[s] == hemisphere[t])]
    report.check(name, "intra-lobe block to elsewhere mean weight ratio",
                 float(np.mean(inside) / np.mean(outside)), 3.921)

    laplacian = np.zeros((size, size))
    for s, t, w in zip(source, target, weight):
        laplacian[s, s] += w
        laplacian[t, t] += w
        laplacian[s, t] -= w
        laplacian[t, s] -= w
    eigenvalues = np.linalg.eigvalsh(laplacian)
    report.check(name, "Fiedler value of L = D - A", float(eigenvalues[1]),
                 0.772254)

    for cells in (1, 2, 4, 8):
        lines = Path(f"data/connectome/fornari83/graph_fem_{cells}.txt") \
            .read_text().split("\n")
        header = lines[0].split()
        rows = [line.split() for line in lines[1:] if line.strip()]
        same = (int(header[0]) == size and int(header[1]) == len(edges)
                and all(int(r[0]) == s and int(r[1]) == t
                        for r, s, t in zip(rows, source, target))
                and all(float(r[2]) == 1.0 and int(r[3]) == cells
                        for r in rows))
        report.check(name, f"graph_fem_{cells}: edge order of edges.csv, "
                           f"unit length, {cells} cells",
                     1.0 if same else 0.0, 1.0)


def check_19_scheme(report):
    """Lobe separation at unit scale under every stored discretization."""
    name = "19 scheme_sensitivity"
    rows = read_csv(BENCH / "19_fisher_kolmogorov_fornari83/results"
                            "/scheme_sensitivity.csv")
    spreads = {(row["model"], row["scheme"], float(row["dt"]),
                int(row["cells_per_edge"])): float(row["lobe_spread_years"])
               for row in rows}
    report.check(name, "rows", float(len(rows)), 9.0)
    report.check(name, "nodal spread", spreads[("nodal", "backward_euler",
                                                0.4, 1)], 3.8703e-07)
    report.check(name, "FEM semi-implicit spread, one cell",
                 spreads[("fem", "corti_semi_implicit", 0.4, 1)], 7.576e-03)
    report.check(name, "FEM backward Euler spread, dt 0.05, eight cells",
                 spreads[("fem", "backward_euler", 0.05, 8)], 3.848e-02)
    report.check(name, "largest spread over every discretization",
                 max(spreads.values()), 3.902e-02)
    report.check(name, "every spread below 0.04 years",
                 1.0 if max(spreads.values()) < 0.04 else 0.0, 1.0)
    # The chapter states the range of the absolute network crossing over
    # the finite element runs of the table, 11.50 to 13.80 years.
    crossings = [float(row["t50_network_years"]) for row in rows
                 if row["model"] == "fem"]
    report.check(name, "earliest network crossing over the runs",
                 min(crossings), 11.5044, "years")
    report.check(name, "latest network crossing over the runs",
                 max(crossings), 13.8032, "years")


def check_24_views(report):
    """Data behind the four-view brain-network figure, after Fornari fig. 5.

    The figure colours and sizes the connections by their fibre number on a
    linear ramp between the smallest and largest value; these checks pin the
    endpoints of that ramp, the counts, and the reproduction of the published
    fibre statistics at printed precision (1 to 596, mean 40.2).
    """
    name = "24 connectome_views"
    edges = read_csv(Path("data/connectome/fornari83/edges.csv"))
    fibres = [float(row["fibre_number"]) for row in edges]
    lengths = [float(row["fibre_length_mm"]) for row in edges]
    report.check(name, "connections drawn", float(len(edges)), 1130.0)
    report.check(name, "smallest fibre number", min(fibres), 1.0)
    report.check(name, "largest fibre number", max(fibres), 595.5)
    report.check(name, "mean fibre number", sum(fibres) / len(fibres),
                 40.1619)
    report.check(name, "published mean fibre number, printed precision",
                 round(sum(fibres) / len(fibres), 1), 40.2)
    report.check(name, "mean fibre length", sum(lengths) / len(lengths),
                 38.4009, "mm")
    report.check(name, "published mean fibre length, printed precision",
                 round(sum(lengths) / len(lengths), 2), 38.40, "mm")


def check_18_timestep(report):
    """Numbers the time-step study reports.

    The front positions and the sentinel are shown on the figure; the two-point
    slopes are recorded in the benchmark README only, deliberately not drawn,
    and are checked here so that the README cannot drift from the data.
    """
    name = "18 time_step_study"
    rows = read_csv(BENCH / "18_fisher_kolmogorov_1d_sensitivity/results"
                           "/time_step_study.csv")
    table = {float(row["dt"]): row for row in rows}
    reference = float(table[0.025]["front_position"])
    for label, key, divisor, expected in (
            ("L2", "l2_error", 1.15408, 1.438),
            ("Linf", "max_error", 1.0, 1.048)):
        coarse = float(table[0.1][key]) / divisor
        fine = float(table[0.05][key]) / divisor
        report.check(name, f"two-point slope, {label}",
                     math.log(coarse / fine) / math.log(2.0), expected)
    coarse = abs(float(table[0.1]["front_position"]) - reference) / reference
    fine = abs(float(table[0.05]["front_position"]) - reference) / reference
    report.check(name, "two-point slope, front position",
                 math.log(coarse / fine) / math.log(2.0), 1.728)
    # The sentinel that must not be read as a measured front position.
    for step in (0.2, 0.3, 0.4):
        report.check(name, f"front sentinel at dt={step:g}",
                     float(table[step]["front_position"]), 1.0)

    # The quantities the figure actually draws, recomputed from the profiles.
    import numpy as np
    profiles = {}
    for row in read_csv(BENCH / "18_fisher_kolmogorov_1d_sensitivity/results"
                                "/time_step_profiles.csv"):
        profiles.setdefault(float(row["dt"]), []).append(
            (float(row["x"]), float(row["c"])))
    profiles = {step: np.array(sorted(points))
                for step, points in profiles.items()}
    reference_profile = profiles[0.025]

    def crossings(profile, level=0.5):
        x, c = profile[:, 0], profile[:, 1]
        return [x[k - 1] + (level - c[k - 1]) / (c[k] - c[k - 1])
                * (x[k] - x[k - 1])
                for k in range(1, len(c))
                if (c[k - 1] - level) * (c[k] - level) < 0]

    with_front = [step for step in profiles if crossings(profiles[step])]
    report.check(name, "runs that cross c = 0.5", float(len(with_front)), 3.0)
    for step in (0.2, 0.3, 0.4):
        report.check(name, f"no crossing at dt={step:g}",
                     float(len(crossings(profiles[step]))), 0.0)
    report.check(name, "minimum concentration at dt=0.2",
                 float(profiles[0.2][:, 1].min()), 0.9320)

    for step, expected in ((0.05, 0.106268), (0.1, 0.288002),
                           (0.4, 0.542938)):
        difference = profiles[step][:, 1] - reference_profile[:, 1]
        rms = float(np.sqrt(np.trapezoid(difference ** 2,
                                         profiles[step][:, 0]) / 2.0))
        report.check(name, f"e_2 at dt={step:g}", rms, expected)
        # e_2 must be exactly the stored L2 column divided by sqrt(2).
        report.check(name, f"e_2 equals stored L2 / sqrt(2) at dt={step:g}",
                     rms, float(table[step]["l2_error"]) / math.sqrt(2.0))


def check_19(report):
    """Lobe separation and refinement rates of the Fornari comparison."""
    name = "19 biomarker_comparison"
    base = BENCH / "19_fisher_kolmogorov_fornari83/results"
    spreads = {}
    networks = {}
    for label, filename, expected in (
            ("nodal", "nodal_biomarkers.csv", 3.8703e-07),
            ("FEM", "fem_biomarkers.csv", 8.378e-03)):
        rows = read_csv(base / filename)
        times = [float(row["time"]) for row in rows]
        crossings = [crossing(times, [float(row[lobe]) for row in rows], 50.0)
                     for lobe in ("temporal", "frontal", "parietal",
                                  "occipital")]
        spreads[label] = max(crossings) - min(crossings)
        report.check(name, f"{label} lobe spread", spreads[label], expected,
                     "years")
        networks[label] = crossing(
            times, [float(row["global"]) for row in rows], 50.0)

    # The comparison prose of the chapter states both spreads. Derive the
    # two-significant-digit LaTeX strings from the stored curves and require
    # them verbatim, so the sentence cannot drift from the stored run, as it
    # once did after the weight change.
    chapter = Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8")
    for label, spread in spreads.items():
        exponent = math.floor(math.log10(spread))
        mantissa = round(spread / 10.0 ** exponent, 1)
        if mantissa >= 10.0:
            mantissa, exponent = mantissa / 10.0, exponent + 1
        needle = f"{mantissa:.1f}\\cdot10^{{{exponent}}}"
        report.check_contains(name, f"chapter states the {label} spread",
                              chapter, needle)
    for label, value in networks.items():
        report.check_contains(name,
                              f"chapter states the {label} network crossing",
                              chapter, f"{value:.2f}")

    # The reproduction-limits section states that the coincidence at unit
    # scale does not depend on how the regions are grouped: all 83 regional
    # curves cross within a stated spread. Recompute it from the stored
    # per-vertex profiles.
    for label, filename, expected in (
            ("nodal", "nodal_profiles.csv", 6.404e-05),
            ("FEM", "fem_profiles.csv", 0.1998)):
        rows = read_csv(base / filename)
        times = [float(row["time"]) for row in rows]
        keys = [key for key in rows[0] if key.startswith("node_")]
        vertex = [crossing(times, [100.0 * float(row[key]) for row in rows],
                           50.0) for key in keys]
        report.check(name, f"{label} 83-vertex crossing spread",
                     max(vertex) - min(vertex), expected, "years")
    report.check_contains(name, "chapter states the vertex-level spreads",
                          chapter, "6.4\\cdot10^{-5}")
    report.check_contains(name, "chapter states the FEM vertex spread",
                          chapter, "$0.20$")

    name = "19 refinement"
    rows = read_csv(base / "space_refinement.csv")
    for row, expected in zip(rows, (0.4992, 0.1287, 0.0246, 0.0)):
        report.check(name, f"max difference, {row['cells_per_edge']} cells",
                     float(row["max_biomarker_difference"]), expected, "%")
    rows = read_csv(base / "time_refinement.csv")
    rates = [float(row["rate"]) for row in rows if row["rate"]]
    for rate, expected in zip(rates, (1.0556, 1.0262, 1.0127)):
        report.check(name, "observed temporal rate", rate, expected)


def check_20(report):
    """Crossing times of the conversion-rate study."""
    name = "20 alpha_sensitivity"
    rows = read_csv(BENCH / "20_fisher_kolmogorov_alpha_sensitivity/results"
                           "/alpha_sensitivity.csv")
    expected = {("Nodal reference", 0.5): 11.1060671084,
                ("Nodal reference", 0.1): 59.3308235145,
                ("Metric-graph FEM", 0.5): 12.6792077931,
                ("Metric-graph FEM", 0.1): 67.8076795244}
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
    expected = {"frontal": 0.202253, "temporal": 0.167108,
                "parietal": 0.096674, "insular": 0.129048,
                "limbic": 0.210617, "occipital": 0.072868,
                "subcortical": 0.127734}
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

    # The report states that the ordering does not rest on the chosen
    # normalization: over six weight scales the agreement stays between 18
    # and 21 of 21, occipital last and frontal/limbic first at every scale.
    name = "21 weight_scale"
    rows = read_csv(BENCH / "21_fisher_kolmogorov_corti83/results"
                            "/weight_scale_sensitivity.csv")
    report.check(name, "scales tested", float(len(rows)), 6.0)
    agreements = [int(row["agreeing_pairs_of_21"]) for row in rows]
    report.check(name, "smallest agreement", float(min(agreements)), 18.0)
    report.check(name, "largest agreement", float(max(agreements)), 21.0)
    damkohler = [float(row["damkohler"]) for row in rows]
    report.check(name, "smallest Da (literal scale)", min(damkohler), 0.1621)
    report.check(name, "largest Da", max(damkohler), 57.27)
    rankings = [row["ranking_at_final_time"].split(" > ") for row in rows]
    report.check(name, "occipital last at every scale",
                 1.0 if all(r[-1] == "occipital" for r in rankings) else 0.0,
                 1.0)
    report.check(name, "frontal and limbic first two at every scale",
                 1.0 if all(set(r[:2]) == {"frontal", "limbic"}
                            for r in rankings) else 0.0, 1.0)
    unit = next(row for row in rows if row["weight_scale"] == "1")
    report.check(name, "scale 1 agreement equals the stored run",
                 float(unit["agreeing_pairs_of_21"]), float(concordant))


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
    report.check(name, "Fiedler value", summary["fiedler_value"], 0.772254)
    report.check(name, "benchmark 21 scaling", summary["benchmark_21_scaling"],
                 1.0 / summary["maximum_adjacency"])
    rows = read_csv(base / "diffusion_scaling.csv")
    for row in rows:
        scaling = float(row["diffusion_scaling"])
        report.check(name, f"Da at rho={scaling:g}", float(row["damkohler"]),
                     summary["alpha"] / (scaling * summary["fiedler_value"]))
    bounded = [float(row["diffusion_scaling"]) for row in rows
               if row["fem_corti_crank_nicolson"] == "bounded"
               and row["fem_backward_euler"] == "bounded"]
    unbounded = [float(row["diffusion_scaling"]) for row in rows
                 if row["fem_corti_crank_nicolson"] == "unbounded"
                 and row["fem_backward_euler"] == "unbounded"]
    report.check(name, "sweep rows", float(len(rows)), 13.0)
    # The validity boundary the report states: bounded up to rho = 0.05
    # (Da = 12.9), unbounded from rho = 0.04 (Da = 16.2) onwards, with the
    # two schemes agreeing at every scaling.
    report.check(name, "smallest bounded scaling", min(bounded), 0.05)
    report.check(name, "largest unbounded scaling", max(unbounded), 0.04)
    report.check(name, "every row classified alike by both schemes",
                 float(len(bounded) + len(unbounded)), float(len(rows)))
    damkohler = {float(row["diffusion_scaling"]): float(row["damkohler"])
                 for row in rows}
    report.check(name, "Da of the last bounded scaling", damkohler[0.05],
                 12.949)
    report.check(name, "Da of the first unbounded scaling", damkohler[0.04],
                 16.186)
    spread = {float(row["diffusion_scaling"]): float(row["lobe_spread_years"])
              for row in rows}
    report.check(name, "lobe spread at rho=1", spread[1.0], 3.8703e-07,
                 "years")
    report.check(name, "lobe spread at rho=0.001", spread[0.001], 6.9863,
                 "years")

    # The stabilized variant the report describes: with the lumped mass the
    # finite element solution stays within [0,1] at every scaling, with both
    # schemes and at four elements per connection where the consistent mass
    # diverged; at the literal scale it agrees with the consistent mass.
    name = "23 fem_lumped"
    lumped = read_csv(base / "fem_lumped_sweep.csv")
    report.check(name, "runs", float(len(lumped)), 28.0)
    report.check(name, "smallest transient minimum",
                 min(float(row["transient_min"]) for row in lumped), 0.0)
    report.check(name, "largest transient maximum",
                 max(float(row["transient_max"]) for row in lumped), 1.0)
    report.check(name, "every run within [0,1]",
                 1.0 if all(-1e-12 <= float(row["transient_min"])
                            and float(row["transient_max"]) <= 1.0 + 1e-12
                            for row in lumped) else 0.0, 1.0)
    unit = next(row for row in lumped if row["diffusion_scaling"] == "1.0"
                and row["scheme"] == "be_lumped")
    report.check(name, "lumped BE network crossing at rho=1",
                 float(unit["t50_network_years"]), 12.7215, "years")
    report.check(name, "lumped BE lobe spread at rho=1",
                 float(unit["lobe_spread_years"]), 2.679e-02, "years")
    refined = [row for row in lumped if row["cells_per_edge"] == "4"]
    report.check(name, "four-element runs at rho=0.005", float(len(refined)),
                 2.0)

    # The transient extremes the validity-boundary discussion quotes: the
    # violation grows with Da before the failure, and at Da = 12.9 the
    # solution leaves [0,1] and still recovers.
    rows = read_csv(BENCH / "19_fisher_kolmogorov_fornari83/results"
                           "/fem_biomarkers.csv")
    report.check(name, "transient minimum at Da = 0.65",
                 min(float(row["min"]) for row in rows), -1.337e-4)
    rows = read_csv(base / "fem_transient_rho_0p05.csv")
    report.check(name, "transient minimum at Da = 12.9",
                 min(float(row["min"]) for row in rows), -0.16743)
    report.check(name, "transient maximum at Da = 12.9",
                 max(float(row["max"]) for row in rows), 1.04356)
    report.check(name, "recovers to the plateau at Da = 12.9",
                 float(rows[-1]["min"]), 0.999424)


def check_23_stabilization(report):
    """The three-row separation figure: nodal, consistent and lumped FEM."""
    name = "23 diffusion_scaling rows"
    rows = read_csv(BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
                            "/diffusion_scaling_summary_rows.csv")
    table = {(row["scheme"], row["diffusion_scaling"]): row for row in rows}
    report.check(name, "runs", float(len(rows)), 9.0)
    for scaling, expected in (("1.0", 3.870e-07), ("0.05", 0.6529),
                              ("0.005", 4.380)):
        report.check(name, f"nodal lobe spread at rho={scaling}",
                     float(table[("nodal", scaling)]["lobe_spread_years"]),
                     expected, "years")
        report.check(name, f"nodal within [0,1] at rho={scaling}",
                     1.0 if float(table[("nodal", scaling)]["min_concentration"]) >= -1e-12
                     and float(table[("nodal", scaling)]["max_concentration"]) <= 1 + 1e-12
                     else 0.0, 1.0)
    report.check(name, "consistent, rho=1: minimum",
                 float(table[("be", "1.0")]["min_concentration"]), -1.337e-4)
    report.check(name, "consistent, rho=0.05: minimum",
                 float(table[("be", "0.05")]["min_concentration"]), -0.0798)
    report.check(name, "consistent, rho=0.05: maximum",
                 float(table[("be", "0.05")]["max_concentration"]), 1.0206)
    report.check(name, "consistent, rho=0.05: lobe spread",
                 float(table[("be", "0.05")]["lobe_spread_years"]), 6.18,
                 "years")
    report.check(name, "consistent, rho=0.005: minimum",
                 float(table[("be", "0.005")]["min_concentration"]), -0.4139)
    report.check(name, "consistent, rho=0.005: maximum",
                 float(table[("be", "0.005")]["max_concentration"]), 1.1712)
    report.check(name, "consistent, rho=0.005: run stops at",
                 float(table[("be", "0.005")]["last_stored_time"]), 15.6,
                 "years")
    for scaling, spread in (("1.0", 0.0268), ("0.05", 5.09), ("0.005", 9.37)):
        row = table[("be_lumped", scaling)]
        report.check(name, f"lumped, rho={scaling}: within [0,1]",
                     1.0 if float(row["min_concentration"]) >= -1e-12
                     and float(row["max_concentration"]) <= 1.0 + 1e-12
                     else 0.0, 1.0)
        report.check(name, f"lumped, rho={scaling}: reaches T=40",
                     float(row["last_stored_time"]), 40.0, "years")
        report.check(name, f"lumped, rho={scaling}: lobe spread",
                     float(row["lobe_spread_years"]), spread, "years")

    # The Damkohler and boundary prose of the chapter quote the FEM spreads
    # at rho = 0.05 (consistent) and 0.005 (lumped) to one decimal; derive
    # the strings from the summary and require them verbatim.
    chapter = Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8")
    for scheme, scaling, what in (("be", "0.05", "consistent"),
                                  ("be_lumped", "0.005", "lumped")):
        value = float(table[(scheme, scaling)]["lobe_spread_years"])
        report.check_contains(name, f"chapter states the {what} FEM spread",
                              chapter, f"${value:.1f}$ years")


def check_23_mass_spectrum(report):
    """The Damkohler number of the nodal model against the finite element one:
    the Fiedler value of L against the identity, the lumped and the consistent
    mass matrix, recomputed from the stored edge list."""
    import numpy as np
    name = "23 mass spectrum"
    edges = read_csv(Path("data/connectome/fornari83/edges.csv"))
    n = 1 + max(max(int(e["source"]), int(e["target"])) for e in edges)
    laplacian = np.zeros((n, n))
    lumped = np.zeros(n)
    consistent = np.zeros((n, n))
    for edge in edges:
        i, j = int(edge["source"]), int(edge["target"])
        w = float(edge["connectivity_weight"])
        laplacian[i, i] += w
        laplacian[j, j] += w
        laplacian[i, j] -= w
        laplacian[j, i] -= w
        lumped[i] += 0.5
        lumped[j] += 0.5
        consistent[i, i] += 1.0 / 3.0
        consistent[j, j] += 1.0 / 3.0
        consistent[i, j] += 1.0 / 6.0
        consistent[j, i] += 1.0 / 6.0

    def fiedler(mass=None):
        reduced = laplacian
        if mass is not None:
            inverse = np.linalg.inv(np.linalg.cholesky(mass))
            reduced = inverse @ laplacian @ inverse.T
        return float(np.sort(np.linalg.eigvalsh(reduced))[1])

    recomputed = {"identity": fiedler(), "lumped": fiedler(np.diag(lumped)),
                  "consistent": fiedler(consistent)}
    rows = {row["mass_matrix"]: row for row in read_csv(
        BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
                "/mass_spectrum.csv")}
    for mass, expected in (("identity", 0.772254), ("lumped", 0.0606268),
                           ("consistent", 0.0635326)):
        report.check(name, f"lambda_2 with the {mass} mass, stored",
                     float(rows[mass]["lambda_2"]), expected)
        report.check(name, f"lambda_2 with the {mass} mass, recomputed",
                     recomputed[mass], expected)
    report.check(name, "lumped diagonal: mean", float(lumped.mean()), 13.6145)
    report.check(name, "lumped diagonal: min", float(lumped.min()), 3.0)
    report.check(name, "lumped diagonal: max", float(lumped.max()), 24.0)
    names = {int(row["node_id"]): row["name"] for row in read_csv(
        Path("data/connectome/fornari83/nodes.csv"))}
    degree = {i: 2.0 * lumped[i] for i in range(n)}
    least = min(degree, key=degree.get)
    report.check(name, "least connected vertex is a frontal pole",
                 1.0 if "frontalpole" in names[least] else 0.0, 1.0)
    caudate = max(degree[i] for i in names if "Caudate" in names[i])
    report.check(name, "a caudate has the largest degree", caudate,
                 max(degree.values()))
    ratio_lumped = recomputed["identity"] / recomputed["lumped"]
    ratio_consistent = recomputed["identity"] / recomputed["consistent"]
    report.check(name, "ratio of Fiedler values, lumped", ratio_lumped, 12.74)
    report.check(name, "ratio of Fiedler values, consistent",
                 ratio_consistent, 12.16)

    # The Damkohler paragraph of the chapter states these numbers; derive
    # its strings from the recomputed values and the stored sweep.
    chapter = Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8")
    words = {11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen"}
    for label, needle in (
            ("lumped lambda_2", f"{recomputed['lumped']:.4f}$"),
            ("consistent lambda_2", f"{recomputed['consistent']:.4f}$"),
            ("identity lambda_2", f"{recomputed['identity']:.4f}$"),
            ("mean lumped mass", f"${lumped.mean():.1f}$ on average"),
            ("lumped mass of the frontal pole",
             f"${lumped.min():g}$ for the frontal pole, which has "
             f"${2 * lumped.min():g}$ connections"),
            ("lumped mass of the caudate",
             f"${lumped.max():g}$ for the caudate, which has "
             f"${2 * lumped.max():g}$"),
            ("ratio range", f"{words[math.floor(ratio_consistent)]} to "
                            f"{words[math.ceil(ratio_lumped)]} times"),
            ("FEM Da at rho = 1",
             f"$\\mathrm{{Da}}={0.5 / recomputed['consistent']:.1f}$")):
        report.check_contains(name, f"chapter states the {label}", chapter,
                              needle)
    # The Damkohler number of the finite element model at rho = 0.05 and
    # 0.005, and the nodal spreads that bracket the former in the stored sweep.
    effective = 0.5 / (0.05 * recomputed["consistent"])
    effective_weak = 0.5 / (0.005 * recomputed["consistent"])
    report.check_contains(name, "chapter states the FEM Da at rho = 0.05 "
                          "and 0.005", chapter,
                          f"\\mathrm{{Da}} \\approx {round(effective, -1):.0f}$ "
                          f"and ${round(effective_weak, -2):.0f}$")
    sweep = sorted(((float(r["damkohler"]), float(r["lobe_spread_years"]))
                    for r in read_csv(
                        BENCH / "23_fisher_kolmogorov_diffusion_scaling"
                                "/results/diffusion_scaling.csv")))
    below = max(s for d, s in sweep if d < effective)
    above = min(s for d, s in sweep if d > effective)
    report.check_contains(name, "chapter states the bracketing nodal spreads",
                          chapter, f"between ${below:.1f}$ and ${above:.1f}$ years")
    report.check(name, "chapter never says Fiedler",
                 float(chapter.count("Fiedler")), 0.0)


def check_27(report):
    """Stage times and expected orderings of the two clinical seedings."""
    name = "27 seeding_patterns"
    base = BENCH / "27_connectome_seeding_patterns/results"
    import numpy as np
    group = {int(row["node_id"]): row["region"]
             for row in read_csv(BENCH / "21_fisher_kolmogorov_corti83"
                                         "/results/reaction_coefficients.csv")}
    # Lumped-mass finite element runs at rho = 0.005 (benchmark 27).
    expected_stages = {"tau": (15.2, 21.6, 27.6),
                       "amyloid": (0.8, 4.8, 11.2)}
    for case in ("tau", "amyloid"):
        rows = read_csv(base / f"{case}_profiles.csv")
        times = np.array([float(row["time"]) for row in rows])
        values = np.array([[float(row[f"node_{k}"]) for k in range(83)]
                           for row in rows])
        report.check(name, f"{case} stays in [0,1]",
                     1.0 if values.min() >= 0.0
                     and values.max() <= 1.0 + 1e-10 else 0.0, 1.0)
        means = values.mean(axis=1)
        for target, expected in zip((0.10, 0.40, 0.80),
                                    expected_stages[case]):
            stage = times[int(np.argmax(means >= target))]
            report.check(name, f"{case} stage at mean {target:.0%}",
                         float(stage), expected, "years")

    # The amyloid ordering the report states: seeded lobes first, then the
    # insula, the limbic belt, and the subcortical nuclei last.
    rows = read_csv(base / "amyloid_profiles.csv")
    times = [float(row["time"]) for row in rows]
    activation = {}
    for k in range(83):
        series = [float(row[f"node_{k}"]) for row in rows]
        first = next((times[i] for i, v in enumerate(series) if v >= 0.5),
                     None)
        activation.setdefault(group[k], []).append(first)
    mean_activation = {g: sum(v) / len(v) for g, v in activation.items()}
    report.check(name, "amyloid subcortical mean activation",
                 mean_activation["subcortical"], 12.98, "years")
    report.check(name, "amyloid limbic mean activation",
                 mean_activation["limbic"], 9.78, "years")
    report.check(name, "amyloid insular mean activation",
                 mean_activation["insular"], 9.6, "years")
    latest_cortical = max(mean_activation[g] for g in
                          ("frontal", "temporal", "parietal", "occipital"))
    report.check(name, "subcortical last, after every cortical lobe",
                 1.0 if mean_activation["subcortical"] > latest_cortical
                 and mean_activation["limbic"] > latest_cortical else 0.0,
                 1.0)

    # The tau row starts at the entorhinal seeds: the first two vertices to
    # reach c = 0.5 must be the entorhinal cortices.
    rows = read_csv(base / "tau_profiles.csv")
    times = [float(row["time"]) for row in rows]
    nodes = {int(row["node_id"]): row["name"]
             for row in read_csv(Path("data/connectome/fornari83/nodes.csv"))}
    firsts = []
    for k in range(83):
        series = [float(row[f"node_{k}"]) for row in rows]
        first = next((times[i] for i, v in enumerate(series) if v >= 0.5),
                     float("inf"))
        firsts.append((first, nodes[k]))
    earliest = sorted(firsts)[:2]
    report.note(name, "first two activated vertices, tau",
                ", ".join(name for _, name in earliest)
                + " (expected the two entorhinal cortices)")
    if not all("entorhinal" in name for _, name in earliest):
        report.failed += 1

    # The staging renders and the diffusion-scaling curves come from the
    # same model at the same scaling, so the tau lobe order behind the
    # staged row must be the one benchmark 23 reports: temporal first,
    # then occipital, parietal, and the weakly connected frontal lobe
    # last. This pins the two figures to each other.
    rows = read_csv(base / "tau_biomarkers.csv")
    times = [float(row["time"]) for row in rows]
    lobe_crossings = {
        lobe: crossing(times, [float(row[lobe]) for row in rows], 50.0)
        for lobe in ("temporal", "frontal", "parietal", "occipital")}
    for lobe, expected in (("temporal", 16.50), ("occipital", 21.65),
                           ("parietal", 22.95), ("frontal", 25.87)):
        report.check(name, f"tau {lobe} lobe crossing at rho=0.005",
                     lobe_crossings[lobe], expected, "years")
    ordered = sorted(lobe_crossings, key=lobe_crossings.get)
    report.check(name, "tau lobe order matches the diffusion scaling",
                 1.0 if ordered == ["temporal", "occipital", "parietal",
                                    "frontal"] else 0.0, 1.0)


def check_19_accuracy(report):
    """The formulation gap the convergence section quantifies."""
    name = "19 accuracy"
    base = BENCH / "19_fisher_kolmogorov_fornari83/results"
    rows = read_csv(base / "nodal_biomarkers.csv")
    times = [float(row["time"]) for row in rows]
    nodal = crossing(times, [float(row["global"]) for row in rows], 50.0)
    rows = read_csv(base / "space_refinement.csv")
    fem = {int(row["cells_per_edge"]): float(row["global_t50"])
           for row in rows}
    report.check(name, "nodal crossing behind the resolved FEM",
                 fem[8] - nodal, 1.605, "years")
    report.check(name, "coarsest FEM behind the resolved FEM",
                 fem[8] - fem[1], 0.0296, "years")

    # The metric-support effect the comparison section describes: the same
    # seed is 0.24% of the nodal average and 0.10% of the metric length of
    # the finite element domain, while the vertex averages of the two
    # models start equal; both finite element averages cross 50% together.
    name = "19 metric_support"
    nodal_rows = read_csv(base / "nodal_biomarkers.csv")
    fem_rows = read_csv(base / "fem_biomarkers.csv")
    metric_rows = read_csv(base / "fem_metric_mass.csv")
    report.check(name, "nodal vertex average at t=0",
                 float(nodal_rows[0]["global"]), 0.2410, "%")
    report.check(name, "FEM vertex average at t=0",
                 float(fem_rows[0]["global"]), 0.2410, "%")
    report.check(name, "FEM metric average at t=0",
                 float(metric_rows[0]["metric_global"]), 0.1018, "%")
    report.check(name, "FEM metric average at t=0 equals seed hats "
                       "over total length",
                 float(metric_rows[0]["metric_global"]),
                 100.0 * 0.1 * (10 + 13) / 2.0 / 1130.0, "%")
    metric_times = [float(row["time"]) for row in metric_rows]
    report.check(name, "FEM metric average crossing",
                 crossing(metric_times,
                          [float(row["metric_global"]) for row in metric_rows],
                          50.0), 12.6792, "years")


def check_25(report):
    """Extremes and spans of the seeding study."""
    name = "25 seeding_vulnerability"
    base = BENCH / "25_connectome_seeding_vulnerability/results"
    # Lumped-mass finite element runs (fully implicit), see benchmark 25.
    for scaling, span, fastest, slowest in (
            (1, 3.9136, "Right-Caudate", "ctx-rh-frontalpole"),
            (0.02, 11.180, "Right-Caudate", "ctx-rh-temporalpole")):
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


def check_26_order(report):
    """Group activation of the calibrated and uniform Corti variants."""
    name = "26 activation_order"
    rows = read_csv(BENCH / "26_connectome_progression/results"
                            "/activation_order_groups.csv")
    means = {(row["variant"], row["region"]):
             float(row["mean_activation_years"]) for row in rows}
    for (variant, region), expected in (
            (("regional", "frontal"), 30.15),
            (("regional", "occipital"), 43.00),
            (("uniform", "occipital"), 35.00),
            (("uniform", "frontal"), 36.20)):
        report.check(name, f"{variant} {region} mean activation",
                     means[(variant, region)], expected, "years")

    groups = ("frontal", "temporal", "parietal", "insular", "limbic",
              "occipital", "subcortical")
    regional = {group: means[("regional", group)] for group in groups}
    uniform = {group: means[("uniform", group)] for group in groups}
    report.check(name, "occipital last with the regional rates",
                 1.0 if max(regional, key=regional.get) == "occipital"
                 else 0.0, 1.0)
    report.check(name, "frontal first with the regional rates",
                 1.0 if min(regional, key=regional.get) == "frontal"
                 else 0.0, 1.0)
    report.check(name, "frontal last with the uniform rate",
                 1.0 if max(uniform, key=uniform.get) == "frontal"
                 else 0.0, 1.0)
    report.check(name, "regional spread",
                 max(regional.values()) - min(regional.values()), 12.85,
                 "years")
    report.check(name, "uniform spread",
                 max(uniform.values()) - min(uniform.values()), 2.95,
                 "years")


def main():
    report = Report()
    for check in (check_18, check_18_timestep, check_19, check_19_accuracy,
                  check_19_topology, check_connectome_consistency,
                  check_19_scheme, check_20, check_21, check_22, check_23,
                  check_23_stabilization, check_23_mass_spectrum,
                  check_24_views, check_25, check_26,
                  check_26_order, check_27):
        try:
            check(report)
        except FileNotFoundError as error:
            report.note(check.__name__, "input", f"missing: {error.filename}")
            report.failed += 1
    report.show()
    return 1 if report.failed else 0


if __name__ == "__main__":
    sys.exit(main())
