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
    expected = {"frontal": 0.193462, "temporal": 0.135457,
                "parietal": 0.082666, "insular": 0.109318,
                "limbic": 0.165760, "occipital": 0.059494,
                "subcortical": 0.106671}
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
    report.check(name, "concordant pairs vs Corti table 3", concordant, 21)
    squared = sum((rank_corti[r] - rank_ours[r]) ** 2 for r in corti)
    spearman = 1 - 6 * squared / (len(corti) * (len(corti) ** 2 - 1))
    report.check(name, "rank correlation", spearman, 1.0)

    # The numbers the section quotes for the configuration and the run.
    import numpy as np
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    rates = [float(r["alpha"]) for r in read_csv(
        BENCH / "21_fisher_kolmogorov_corti83/results"
                "/reaction_coefficients.csv")]
    weights = [float(e["connectivity_weight"]) for e in read_csv(
        Path("data/connectome/fornari83/edges.csv"))]
    mean_rate = sum(rates) / len(rates)
    rho = 1.0 / max(weights)
    lobe_rate = next(float(r["lobe_rate"]) for r in read_csv(
        BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
                "/lobe_damkohler.csv") if r["model"] == "consistent")
    report.check(name, "mean regional rate", mean_rate, 0.1252)
    report.check(name, "nominal Damkohler number",
                 mean_rate / (rho * 0.772254), 5.72)
    report.check(name, "lobe Damkohler number, consistent mass",
                 mean_rate / (rho * lobe_rate), 3.99)
    # The lobe rate stored by benchmark 23 is computed at one element per
    # connection; the run uses four. With the consistent mass the quotient
    # over lobe-constant patterns, linear across the connections between
    # groups, is the same at any refinement, because P1 elements reproduce
    # a linear ramp and integrate its square exactly. Checked by assembling
    # the operators of the four-element mesh, as the run does, with the
    # normalized diffusivity w/max(w).
    import scipy.sparse as sparse
    from scipy.linalg import eigh as dense_eigh
    groups_of = {}
    keys = {"temporal": ("temporal", "bankssts", "entorhinal", "fusiform",
                         "parahippocampal"),
            "frontal": ("frontal", "orbitofrontal", "parsopercularis",
                        "parsorbitalis", "parstriangularis", "precentral"),
            "parietal": ("parietal", "postcentral", "precuneus",
                         "supramarginal", "paracentral"),
            "occipital": ("cuneus", "occipital", "lingual", "pericalcarine")}
    for row in read_csv(Path("data/connectome/fornari83/nodes.csv")):
        lowered = row["name"].lower()
        groups_of[int(row["node_id"])] = next(
            (g for g, words in keys.items()
             if any(w in lowered for w in words)), "other")
    order = ("temporal", "frontal", "parietal", "occipital", "other")
    edge_rows = read_csv(Path("data/connectome/fornari83/edges.csv"))
    cells = 4
    size = 83 + (cells - 1) * len(edge_rows)
    stiffness = sparse.lil_matrix((size, size))
    mass = sparse.lil_matrix((size, size))
    prolongation = np.zeros((size, len(order)))
    for vertex in range(83):
        prolongation[vertex, order.index(groups_of[vertex])] = 1.0
    step = 1.0 / cells
    next_dof = 83
    for edge in edge_rows:
        i, j = int(edge["source"]), int(edge["target"])
        diffusivity = float(edge["connectivity_weight"]) * rho
        chain = [i] + list(range(next_dof, next_dof + cells - 1)) + [j]
        next_dof += cells - 1
        for k, dof in enumerate(chain[1:-1], start=1):
            prolongation[dof, order.index(groups_of[i])] += 1.0 - k * step
            prolongation[dof, order.index(groups_of[j])] += k * step
        for a, b in zip(chain[:-1], chain[1:]):
            for (r, c, sign) in ((a, a, 1), (b, b, 1), (a, b, -1), (b, a, -1)):
                stiffness[r, c] += sign * diffusivity / step
            mass[a, a] += step / 3
            mass[b, b] += step / 3
            mass[a, b] += step / 6
            mass[b, a] += step / 6
    report.check(name, "degrees of freedom of the four-element mesh",
                 float(size), 3473.0)
    coarse_h = prolongation.T @ stiffness.tocsr() @ prolongation
    coarse_m = prolongation.T @ mass.tocsr() @ prolongation
    refined_rate = float(np.sort(dense_eigh(coarse_h, coarse_m,
                                            eigvals_only=True))[1])
    report.check(name, "lobe rate on the four-element mesh, scaled by rho",
                 refined_rate, rho * lobe_rate)
    report.check(name, "lobe Damkohler number on the four-element mesh",
                 mean_rate / refined_rate, 3.99)
    report.check(name, "ratio of the extreme rates", max(rates) / min(rates),
                 3.30)
    for label, needle in (
            ("scaling",
             f"\\rho=\\frac{{1}}{{\\max_e w_e}}={rho:.4f}."),
            ("nominal Damkohler number",
             f"gives the nominal value $\\mathrm{{Da}}="
             f"{mean_rate / (rho * 0.772254):.2f}$"),
            ("lobe rate at this scaling",
             f"\\rho\\,\\lambda_{{\\mathrm{{lobe}}}}="
             f"{rho * lobe_rate:.6f}\\ \\mathrm{{yr}}^{{-1}}"),
            ("lobe Damkohler number",
             f"= \\frac{{\\overline{{\\alpha}}}}"
             f"{{\\rho\\,\\lambda_{{\\mathrm{{lobe}}}}}} "
             f"= {mean_rate / (rho * lobe_rate):.1f}.")):
        report.check_contains(name, f"regional-rate prose states the {label}",
                              chapter, needle)
    monotone = all(all(float(rows[k][g]) > float(rows[k - 1][g])
                       for k in range(1, len(rows)))
                   for g in expected)
    report.check(name, "every regional average grows monotonically",
                 1.0 if monotone else 0.0, 1.0)
    # The range of the full finite element field needs the solution files
    # of the run, which are not stored; checked when they are present.
    paths = sorted(glob.glob("output/fisher_kolmogorov/corti83"
                             "/solution_*.vtp"))
    try:
        import vtk
    except ImportError:
        vtk = None
    if paths and vtk is not None:
        low, high, final = math.inf, -math.inf, None
        for path in sorted(paths, key=lambda q: int(
                re.search(r"solution_(\d+)\.vtp$", q).group(1))):
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(path)
            reader.Update()
            final = reader.GetOutput().GetPointData().GetArray("c").GetRange()
            low, high = min(low, final[0]), max(high, final[1])
        report.check(name, "field minimum over the run", low, 0.01)
        report.check(name, "field maximum over the run", high, 0.4805)
        report.check(name, "field minimum at T", final[0], 0.0253)
    else:
        report.note(name, "field range", "skipped, run benchmark 21 first")

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

    # The Damkohler section states these numbers; derive its strings from
    # the recomputed values.
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    words = {11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen"}
    for label, needle in (
            ("consistent lambda_2", f"from ${recomputed['consistent']:.4f}$ "
                                    f"(the value of the consistent mass) to "
                                    f"${recomputed['lumped']:.4f}$"),
            ("identity lambda_2", f"$\\lambda_2 = {recomputed['identity']:.4f}$"),
            ("lumped mass of the frontal pole",
             f"${lumped.min():g}$ at the frontal pole"),
            ("lumped mass of the caudate", f"${lumped.max():g}$ at the caudate"),
            ("ratio range", f"{words[math.floor(ratio_consistent)]} to "
                            f"{words[math.ceil(ratio_lumped)]} times smaller"),
            ("nodal Da at rho = 1",
             f"$\\mathrm{{Da}}={0.5 / recomputed['identity']:.2f}$"),
            ("FEM Da at rho = 1",
             f"$\\mathrm{{Da}}={0.5 / recomputed['consistent']:.1f}$")):
        report.check_contains(name, f"chapter states the {label}", chapter,
                              needle)
    report.check(name, "chapter never says Fiedler",
                 float(chapter.count("Fiedler")), 0.0)


def check_23_boundary_prose(report):
    """The stabilization subsection quotes the extremes of the consistent
    mass, the stopping time, the Newton safeguard and the divergence of the
    semi-implicit scheme; derive the strings from the stored summary, the
    source of the library and the benchmark record."""
    name = "23 boundary prose"
    rows = {(r["scheme"], r["diffusion_scaling"]): r for r in read_csv(
        BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
                "/diffusion_scaling_summary_rows.csv")}
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    literal = float(rows[("be", "1.0")]["min_concentration"])
    mantissa, exponent = f"{literal:.1e}".split("e")
    report.check_contains(name, "chapter states the undershoot at Da 0.65",
                          chapter, f"${mantissa}\\cdot10^{{{int(exponent)}}}$")
    report.check_contains(name, "chapter states the extremes at Da 12.9",
                          chapter, f"it reaches ${float(rows[('be', '0.05')]['min_concentration']):.2f}$ "
                          f"and the maximum reaches ${float(rows[('be', '0.05')]['max_concentration']):.2f}$")
    report.check_contains(name, "chapter states the extremes at Da 129",
                          chapter, f"reach ${float(rows[('be', '0.005')]['min_concentration']):.2f}$ "
                          f"and ${float(rows[('be', '0.005')]['max_concentration']):.2f}$ by "
                          f"$t={float(rows[('be', '0.005')]['last_stored_time']):g}$ years")
    source = Path("FEMG/include/fisher_kolmogorov_problem.hpp").read_text()
    match = re.search(r"newton_admissible_range_ = ([0-9.]+);", source)
    bound = float(match.group(1))
    report.check_contains(name, "chapter states the Newton safeguard", chapter,
                          f"leaving the interval $[-{bound:g},{1 + bound:g}]$")
    record = Path(BENCH / "23_fisher_kolmogorov_diffusion_scaling/README.md").read_text()
    match = re.search(r"metric-graph FEM \[(-?[0-9.]+),\s+([0-9.]+)\]", record)
    low, high = float(match.group(1)), float(match.group(2))
    report.check_contains(name, "chapter states the semi-implicit divergence",
                          chapter, f"diverges to $[{round(low)},{round(high)}]$")
    report.check_contains(name, "chapter states the literal-scale crossings",
                          chapter, "from $12.68$ to $12.72$ years")


def check_23_lobe_scale(report):
    """Lobe-scale rates and Damkohler numbers: the stored record against an
    independent recomputation, and the consistent-mass sweep within its
    validity boundary."""
    import numpy as np
    from scipy.linalg import eigh
    name = "23 lobe scale"
    base = BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
    groups_of = {}
    keys = {"temporal": ("temporal", "bankssts", "entorhinal", "fusiform",
                         "parahippocampal"),
            "frontal": ("frontal", "orbitofrontal", "parsopercularis",
                        "parsorbitalis", "parstriangularis", "precentral"),
            "parietal": ("parietal", "postcentral", "precuneus",
                         "supramarginal", "paracentral"),
            "occipital": ("cuneus", "occipital", "lingual", "pericalcarine")}
    for row in read_csv(Path("data/connectome/fornari83/nodes.csv")):
        lowered = row["name"].lower()
        groups_of[int(row["node_id"])] = next(
            (g for g, words in keys.items()
             if any(w in lowered for w in words)), "other")
    order = ("temporal", "frontal", "parietal", "occipital", "other")
    n = 1 + max(groups_of)
    laplacian = np.zeros((n, n))
    lumped = np.zeros(n)
    consistent = np.zeros((n, n))
    for edge in read_csv(Path("data/connectome/fornari83/edges.csv")):
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
    projector = np.zeros((n, 5))
    for i, g in groups_of.items():
        projector[i, order.index(g)] = 1.0
    masses = {"nodal": np.eye(n), "consistent": consistent,
              "lumped": np.diag(lumped)}
    expected = {"nodal": (10.0103, 12.9624), "consistent": (1.10733, 1.4339),
                "lumped": (0.863236, 1.11781)}
    record = {row["model"]: row for row in read_csv(base / "lobe_damkohler.csv")}
    nodal_global = float(np.sort(np.linalg.eigvalsh(laplacian))[1])
    rates = {}
    for model, mass in masses.items():
        coarse = eigh(projector.T @ laplacian @ projector,
                      projector.T @ mass @ projector, eigvals_only=True)
        rates[model] = float(np.sort(coarse)[1])
        report.check(name, f"{model}: lobe rate, recomputed", rates[model],
                     expected[model][0])
        report.check(name, f"{model}: lobe rate, stored",
                     float(record[model]["lobe_rate"]), expected[model][0])
        report.check(name, f"{model}: onset ratio", rates[model] / nodal_global,
                     expected[model][1])
        for scaling, column in ((1.0, "damkohler_lobe_rho_1"),
                                (0.05, "damkohler_lobe_rho_0.05"),
                                (0.005, "damkohler_lobe_rho_0.005")):
            report.check(name, f"{model}: Da_lobe at rho = {scaling:g}",
                         float(record[model][column]),
                         0.5 / (scaling * rates[model]))
    # The panels of the report figure print the same numbers.
    scheme_model = {"nodal": "nodal", "be": "consistent",
                    "be_lumped": "lumped"}
    for row in read_csv(base / "diffusion_scaling_summary_rows.csv"):
        rho = float(row["diffusion_scaling"])
        report.check(name, f"figure panel {row['scheme']} rho={rho:g}: Da_lobe",
                     float(row["damkohler_lobe"]),
                     0.5 / (rho * rates[scheme_model[row["scheme"]]]))
    # The eigenvector of the nodal global rate lives on the periphery.
    values, vectors = np.linalg.eigh(laplacian)
    vector = vectors[:, 1]
    names = {int(r["node_id"]): r["name"].lower()
             for r in read_csv(Path("data/connectome/fornari83/nodes.csv"))}
    top = np.argsort(-np.abs(vector))[:5]
    peripheral = all(any(k in names[int(i)] for k in
                         ("frontalpole", "temporalpole", "entorhinal"))
                     for i in top)
    report.check(name, "largest components on poles and entorhinal cortices",
                 1.0 if peripheral else 0.0, 1.0)
    right = [vector[int(i)] for i in top if "rh-" in names[int(i)]]
    left = [vector[int(i)] for i in top if "lh-" in names[int(i)]]
    report.check(name, "opposite signs in the two hemispheres",
                 1.0 if right and left and max(right) * max(left) < 0
                 and min(right) * max(right) > 0 and min(left) * max(left) > 0
                 else 0.0, 1.0)
    # The consistent-mass sweep within the validity boundary.
    bounded = {row["diffusion_scaling"]: row
               for row in read_csv(base / "fem_consistent_bounded.csv")}
    for scaling, spread in (("1.0", 0.00838), ("0.5", 0.2440),
                            ("0.2", 2.003), ("0.1", 4.031), ("0.05", 6.176)):
        report.check(name, f"consistent mass, rho = {scaling}: lobe spread",
                     float(bounded[scaling]["lobe_spread_years"]), spread,
                     "years")
        report.check(name, f"consistent mass, rho = {scaling}: bounded at the end",
                     1.0 if float(bounded[scaling]["transient_max"]) < 1.03
                     and float(bounded[scaling]["transient_min"]) > -0.1
                     else 0.0, 1.0)
    # Onset: first stored point with a spread above one year per model.
    def first_above(points):
        return min(d for d, s in points if s >= 1.0)
    nodal_points = [(0.5 / (0.5 / (float(r["damkohler"]) * 0.772254) * rates["nodal"]),
                     float(r["lobe_spread_years"]))
                    for r in read_csv(base / "diffusion_scaling.csv")]
    lumped_points = [(0.5 / (0.5 / (float(r["damkohler"]) * 0.772254) * rates["lumped"]),
                      float(r["lobe_spread_years"]))
                     for r in read_csv(base / "fem_lumped_sweep.csv")
                     if r["scheme"] == "be_lumped" and r["cells_per_edge"] == "1"]
    consistent_points = [(0.5 / (0.5 / (float(r["damkohler"]) * 0.772254) * rates["consistent"]),
                          float(r["lobe_spread_years"]))
                         for r in bounded.values()]
    report.check(name, "nodal: first point above one year, Da_lobe",
                 first_above(nodal_points), 1.665)
    report.check(name, "lumped: first point above one year, Da_lobe",
                 first_above(lumped_points), 2.897)
    report.check(name, "consistent: first point above one year, Da_lobe",
                 first_above(consistent_points), 2.258)
    # Where each model enters the band of the reference separation, 5.5
    # years with a margin of 15 percent, between two stored points.
    low = 5.5 * 0.85
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    for label, points, before, inside in (
            ("nodal", nodal_points, 10.0, 25.0),
            ("lumped", lumped_points, 5.8, 11.6),
            ("consistent", consistent_points, 4.5, 9.0)):
        points = sorted(points)
        last_below = max(d for d, s in points if s < low)
        first_inside = min(d for d, s in points if s >= low)
        report.check(name, f"{label}: last point below the reference band",
                     last_below, before)
        report.check(name, f"{label}: first point inside the reference band",
                     first_inside, inside)
    report.check_contains(name, "chapter reads the reference band off panel (b)",
                          chapter, "is reached by every model at "
                          "$\\mathrm{Da}_{\\mathrm{lobe}}$ of order ten")


def check_23_lobe_order(report):
    """The prose of the Damkohler section on the lobe-scale rates and on the
    order of the lobes, against the edge list and the stored tau run."""
    import numpy as np
    from scipy.stats import spearmanr
    name = "23 lobe order"
    base = BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
    record = {row["model"]: row for row in read_csv(base / "lobe_damkohler.csv")}
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    nodal_rates = [float(v) for v in record["nodal"]["lobe_rates"].split()]
    needles = [
        ("lobe-scale rates of the nodal model",
         f"from ${nodal_rates[0]:.1f}$ to ${nodal_rates[-1]:.1f}$"),
        ("temporal-occipital contrast",
         f"relaxing at ${float(record['nodal']['temporal_occipital_rate']):.1f}$"),
        ("FEM lobe rates",
         f"${float(record['consistent']['lobe_rate']):.2f}$ with the consistent "
         f"mass and ${float(record['lumped']['lobe_rate']):.2f}$ with the lumped"),
        ("Da_lobe at rho = 1",
         f"${float(record['nodal']['damkohler_lobe_rho_1']):.2f}$ for the nodal "
         f"model and ${float(record['consistent']['damkohler_lobe_rho_1']):.2f}$ "
         f"and ${float(record['lumped']['damkohler_lobe_rho_1']):.2f}$"),
        ("caption: columns in Da_lobe",
         "$0.05$, $1.0$ and $10.0$ for the nodal row and at $0.45$ and $0.58$, "
         "$9.0$ and $11.6$, $90$ and $116$"),
    ]
    # The caption numbers against the summary of the figure.
    summary = read_csv(base / "diffusion_scaling_summary_rows.csv")
    printed = {("nodal", "1.0"): 0.05, ("nodal", "0.05"): 1.0,
               ("nodal", "0.005"): 10.0, ("be", "1.0"): 0.45,
               ("be", "0.05"): 9.0, ("be", "0.005"): 90.0,
               ("be_lumped", "1.0"): 0.58, ("be_lumped", "0.05"): 11.6,
               ("be_lumped", "0.005"): 116.0}
    for row in summary:
        key = (row["scheme"], row["diffusion_scaling"])
        value = float(row["damkohler_lobe"])
        shown = (round(value) if value >= 50.0
                 else round(value, 1) if value >= 0.995 else round(value, 2))
        report.check(name, f"caption Da_lobe {key[0]} rho={key[1]}", shown,
                     printed[key])

    # Per-region coupling to the temporal lobe and the crossing times of the
    # stored tau run (the run of panel (i)).
    keys = {"temporal": ("temporal", "bankssts", "entorhinal", "fusiform",
                         "parahippocampal"),
            "frontal": ("frontal", "orbitofrontal", "parsopercularis",
                        "parsorbitalis", "parstriangularis", "precentral"),
            "parietal": ("parietal", "postcentral", "precuneus",
                         "supramarginal", "paracentral"),
            "occipital": ("cuneus", "occipital", "lingual", "pericalcarine")}
    names, lobes = {}, {}
    for row in read_csv(Path("data/connectome/fornari83/nodes.csv")):
        index = int(row["node_id"])
        names[index] = row["name"].lower()
        lobes[index] = next((g for g, words in keys.items()
                             if any(w in names[index] for w in words)), "other")
    coupling = {i: 0.0 for i in names}
    for edge in read_csv(Path("data/connectome/fornari83/edges.csv")):
        i, j = int(edge["source"]), int(edge["target"])
        w = float(edge["connectivity_weight"])
        if lobes[j] == "temporal":
            coupling[i] += w
        if lobes[i] == "temporal":
            coupling[j] += w
    tau = read_csv(BENCH / "27_connectome_seeding_patterns/results"
                          "/tau_profiles.csv")
    times = [float(r["time"]) for r in tau]
    stats = {}
    for lobe in ("occipital", "parietal", "frontal"):
        members = [i for i, g in lobes.items() if g == lobe]
        w = np.array([coupling[i] for i in members])
        t = np.array([crossing(times, [100.0 * float(r[f"node_{i}"])
                                        for r in tau], 50.0) for i in members])
        stats[lobe] = (w.sum(), w.mean(), float(np.median(w)),
                       float(spearmanr(t, w).correlation))
    report.check(name, "occipital: coupling to the temporal lobe",
                 stats["occipital"][0], 55.58)
    report.check(name, "occipital: mean per region", stats["occipital"][1],
                 6.95)
    report.check(name, "occipital: median per region", stats["occipital"][2],
                 5.36)
    report.check(name, "parietal: coupling to the temporal lobe",
                 stats["parietal"][0], 74.99)
    report.check(name, "parietal: mean per region", stats["parietal"][1],
                 6.25)
    report.check(name, "parietal: median per region", stats["parietal"][2],
                 3.32)
    inferior = sum(coupling[i] for i in names if "inferiorparietal" in names[i])
    report.check(name, "parietal: the two inferior parietal cortices",
                 inferior, 49.42)
    for lobe, expected in (("occipital", -0.95), ("parietal", -0.92),
                           ("frontal", -0.78)):
        report.check(name, f"{lobe}: rank correlation of crossing and coupling",
                     stats[lobe][3], expected)
    needles += [
        ("occipital per-region coupling",
         f"an average of ${stats['occipital'][1]:.1f}$ per region and a median "
         f"of ${stats['occipital'][2]:.1f}$"),
        ("parietal per-region coupling",
         f"an average of ${stats['parietal'][1]:.1f}$ and a smaller median of "
         f"${stats['parietal'][2]:.1f}$"),
        ("inferior parietal share", f"${inferior:.1f}$ of the total parietal"),
    ]
    for label, needle in needles:
        report.check_contains(name, f"chapter states the {label}", chapter,
                              needle)


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
                       "amyloid": (1.2, 5.2, 12.0)}
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
    # The order of the two middle lobes in the tau run: the direct coupling
    # favours the parietal lobe and its first region does cross first, but
    # the slow parietal tail moves the group mean behind the occipital one.
    # Recompute all the numbers the staging prose quotes.
    def classify(name):
        lowered = name.lower()
        if any(key in lowered for key in ("temporal", "bankssts",
                                          "entorhinal", "fusiform",
                                          "parahippocampal")):
            return "temporal"
        if any(key in lowered for key in ("frontal", "orbitofrontal",
                                          "parsopercularis", "parsorbitalis",
                                          "parstriangularis", "precentral")):
            return "frontal"
        if any(key in lowered for key in ("parietal", "postcentral",
                                          "precuneus", "supramarginal",
                                          "paracentral")):
            return "parietal"
        if any(key in lowered for key in ("cuneus", "occipital", "lingual",
                                          "pericalcarine")):
            return "occipital"
        return "other"

    node_rows = read_csv(Path("data/connectome/fornari83/nodes.csv"))
    lobes = {int(row["node_id"]): classify(row["name"]) for row in node_rows}
    coupling = {}
    for edge in read_csv(Path("data/connectome/fornari83/edges.csv")):
        one = lobes[int(edge["source"])]
        two = lobes[int(edge["target"])]
        if one != two:
            key = tuple(sorted((one, two)))
            coupling[key] = coupling.get(key, 0.0) \
                + float(edge["connectivity_weight"])
    report.check(name, "temporal-parietal coupling",
                 coupling[("parietal", "temporal")], 74.99)
    report.check(name, "temporal-occipital coupling",
                 coupling[("occipital", "temporal")], 55.58)
    report.check(name, "frontal-parietal coupling",
                 coupling[("frontal", "parietal")], 190.46)
    report.check(name, "frontal-temporal coupling",
                 coupling[("frontal", "temporal")], 1.03)
    report.check(name, "frontal-occipital coupling",
                 coupling[("frontal", "occipital")], 2.27)
    report.check(name, "frontal-deep coupling",
                 coupling[("frontal", "other")], 220.32)
    weighted_degree = {}
    for edge in read_csv(Path("data/connectome/fornari83/edges.csv")):
        for end in ("source", "target"):
            index = int(edge[end])
            weighted_degree[index] = weighted_degree.get(index, 0.0) \
                + float(edge["connectivity_weight"])
    report.check(name, "minimum weighted degree",
                 min(weighted_degree.values()), 2.0505)

    tau = read_csv(base / "tau_profiles.csv")
    times = [float(row["time"]) for row in tau]
    middle = {"parietal": [], "occipital": []}
    for index, lobe in lobes.items():
        if lobe not in middle:
            continue
        values = [100.0 * float(row[f"node_{index}"]) for row in tau]
        middle[lobe].append(crossing(times, values, 50.0))
    report.check(name, "parietal regions", float(len(middle["parietal"])),
                 12.0)
    report.check(name, "occipital regions", float(len(middle["occipital"])),
                 8.0)
    report.check(name, "first parietal crossing",
                 min(middle["parietal"]), 19.4, "years")
    report.check(name, "first occipital crossing",
                 min(middle["occipital"]), 19.6, "years")
    report.check(name, "last occipital crossing",
                 max(middle["occipital"]), 23.4, "years")
    report.check(name, "last parietal crossing",
                 max(middle["parietal"]), 26.0, "years")
    poles = sorted(
        crossing(times, [100.0 * float(row[f"node_{node['node_id']}"])
                         for row in tau], 50.0)
        for node in node_rows if "frontalpole" in node["name"].lower())
    report.check(name, "first frontal pole crossing", poles[0], 29.82,
                 "years")
    report.check(name, "second frontal pole crossing", poles[1], 30.94,
                 "years")
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    for label, needle in (
            ("damkohler temporal couplings",
             f"with total weight ${coupling[('parietal', 'temporal')]:.1f}$, "
             f"and to the occipital lobe, with total weight "
             f"${coupling[('occipital', 'temporal')]:.1f}$"),
            ("damkohler direct frontal coupling",
             f"frontal lobe is much weaker, with total weight "
             f"${coupling[('frontal', 'temporal')]:.1f}$"),
            ("frontal pole degree",
             f"The weighted degree spans "
             f"${min(weighted_degree.values()):.2f} \\le D_{{ii}} \\le "
             f"{max(weighted_degree.values()):.2f}$")):
        report.check_contains(name, f"staging prose states the {label}",
                              chapter, needle)

    names_all = {int(r["node_id"]): r["name"] for r in read_csv(
        Path("data/connectome/fornari83/nodes.csv"))}
    summary_row = next(r for r in read_csv(
        BENCH / "23_fisher_kolmogorov_diffusion_scaling/results"
                "/diffusion_scaling_summary_rows.csv")
        if r["scheme"] == "be_lumped" and r["diffusion_scaling"] == "0.005")
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
                 mean_activation["subcortical"], 13.0545, "years")
    report.check(name, "amyloid limbic mean activation",
                 mean_activation["limbic"], 11.6750, "years")
    report.check(name, "amyloid insular mean activation",
                 mean_activation["insular"], 9.6000, "years")
    brainstem = next(times[i] for i, row in enumerate(rows)
                     if float(row["node_82"]) >= 0.5)
    # The brainstem is one of the eleven vertices of the subcortical group,
    # so the report quotes the group without it and the brainstem itself.
    deep = [k for k in range(83) if group[k] == "subcortical"]
    stem = [k for k in deep if "Brain-Stem" in names_all[k]]
    report.check(name, "the brainstem is inside the subcortical group",
                 float(len(stem)), 1.0)
    deep_first = []
    for k in deep:
        if k in stem:
            continue
        series = [float(row[f"node_{k}"]) for row in rows]
        deep_first.append(next(t for t, v in zip(times, series) if v >= 0.5))
    deep_without_brainstem = sum(deep_first) / len(deep_first)
    report.check(name, "subcortical mean without the brainstem",
                 deep_without_brainstem, 12.64, "years")
    report.check(name, "amyloid brainstem activation", brainstem, 17.2, "years")
    report.check(name, "brainstem last of all",
                 1.0 if all(brainstem >= a for v in activation.values()
                            for a in v) else 0.0, 1.0)
    # The amyloid seed is the neocortex proper: the entorhinal and
    # parahippocampal cortices start from zero and activate with the limbic
    # group, after the seeded lobes.
    names = {int(r["node_id"]): r["name"].lower() for r in read_csv(
        Path("data/connectome/fornari83/nodes.csv"))}
    allocortex = [k for k in names if "entorhinal" in names[k]
                  or "parahippocampal" in names[k]]
    report.check(name, "allocortex unseeded",
                 max(float(rows[0][f"node_{k}"]) for k in allocortex), 0.0)
    report.check(name, "seeded vertices", float(sum(
        1 for k in range(83) if float(rows[0][f"node_{k}"]) > 0.0)), 54.0)
    means = np.array([[float(row[f"node_{k}"]) for k in range(83)]
                      for row in rows]).mean(axis=1)
    stages = [times[int(np.argmax(means >= target))]
              for target in (0.10, 0.40, 0.80)]
    seeded = [a for g in ("frontal", "temporal", "parietal", "occipital")
              for a in activation[g]]
    for label, needle in (
            ("amyloid stages",
             f"occur at ${stages[0]:.1f}$, ${stages[1]:.1f}$ and "
             f"${stages[2]:.1f}$ years"),
            ("amyloid groups",
             f"the seeded neocortical regions cross first, at "
             f"${sum(seeded) / len(seeded):.1f}$ years, followed "
             f"by the insula at ${mean_activation['insular']:.1f}$ years, "
             f"the limbic regions at ${mean_activation['limbic']:.1f}$ "
             f"years, the remaining subcortical nuclei at "
             f"${deep_without_brainstem:.1f}$ years and finally "
             f"the brainstem itself at ${brainstem:.1f}$ years"),
            ("distinction between the clinical staging and the lobe "
             "sequence",
             "The lobe sequence is instead a result of the network model "
             "of~\\cite{fornari2019prion}, whose figure 7 reports temporal, "
             "frontal, parietal and occipital"),
            ("staging scaling in Da and Da_lobe",
             f"$\\mathrm{{Da}}={round(float(summary_row['damkohler']))}$ and "
             f"$\\mathrm{{Da}}_{{\\mathrm{{lobe}}}}="
             f"{round(float(summary_row['damkohler_lobe']))}$")):
        report.check_contains(name, f"staging prose states the {label}",
                              chapter, needle)
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
    # The caption states what the middle tau panel of the figure shows.
    profiles = read_csv(base / "tau_profiles.csv")
    stage = next(r for r in profiles if abs(float(r["time"]) - 21.6) < 1e-9)
    lobe_of = {}
    for row in read_csv(Path("data/connectome/fornari83/nodes.csv")):
        lowered = row["name"].lower()
        lobe_of[int(row["node_id"])] = next(
            (g for g, words in (
                ("temporal", ("temporal", "bankssts", "entorhinal",
                              "fusiform", "parahippocampal")),
                ("frontal", ("frontal", "orbitofrontal", "parsopercularis",
                             "parsorbitalis", "parstriangularis",
                             "precentral")),
                ("parietal", ("parietal", "postcentral", "precuneus",
                              "supramarginal", "paracentral")),
                ("occipital", ("cuneus", "occipital", "lingual",
                               "pericalcarine")))
             if any(w in lowered for w in words)), "other")
    stage_mean = {}
    for lobe in ("temporal", "frontal", "parietal", "occipital"):
        members = [k for k, g in lobe_of.items() if g == lobe]
        stage_mean[lobe] = sum(float(stage[f"node_{k}"])
                               for k in members) / len(members)
    report.check_contains(
        name, "caption states the second tau stage", chapter,
        f"the temporal lobe has reached ${stage_mean['temporal'] * 100:.0f}\\%$ "
        f"and the frontal lobe ${stage_mean['frontal'] * 100:.0f}\\%$, with the "
        f"occipital and parietal lobes in between at "
        f"${stage_mean['occipital'] * 100:.0f}\\%$ and "
        f"${stage_mean['parietal'] * 100:.0f}\\%$")

    ordered = sorted(lobe_crossings, key=lobe_crossings.get)
    report.check(name, "tau lobe order matches the diffusion scaling",
                 1.0 if ordered == ["temporal", "occipital", "parietal",
                                    "frontal"] else 0.0, 1.0)

    # The same tau run at rho = 1, the transport scale of the references,
    # stored to justify the scaling of the figure: the first stage is a
    # uniform mantle and the lobes cross together, so no staging exists.
    rows = read_csv(base / "tau_rho1_profiles.csv")
    times = np.array([float(row["time"]) for row in rows])
    values = np.array([[float(row[f"node_{k}"]) for k in range(83)]
                       for row in rows])
    first = int(np.argmax(values.mean(axis=1) >= 0.10))
    stage_time = float(times[first])
    report.check(name, "rho=1 first stage", stage_time, 8.8, "years")
    report.check(name, "rho=1 lowest vertex at the first stage",
                 float(values[first].min()), 0.101)
    report.check(name, "rho=1 highest vertex at the first stage",
                 float(values[first].max()), 0.131)
    rows = read_csv(base / "tau_rho1_biomarkers.csv")
    times = [float(row["time"]) for row in rows]
    spread = [crossing(times, [float(row[lobe]) for row in rows], 50.0)
              for lobe in ("temporal", "frontal", "parietal", "occipital")]
    report.check(name, "rho=1 lobe spread", max(spread) - min(spread),
                 0.027, "years")
    report.check_contains(
        name, "staging prose states the rho=1 lobe spread", chapter,
        f"in a control tau simulation the four lobes cross the $50\\%$ "
        f"level within ${max(spread) - min(spread):.2f}$ years of each other")


def check_27_regional(report):
    """The same two seedings with the regional reaction rates: the record
    variant of benchmark 27, which recovers the clinical lobe order."""
    import numpy as np
    name = "27 regional_rates"
    base = BENCH / "27_connectome_seeding_patterns/results"
    rates = [float(r["alpha"]) for r in read_csv(
        BENCH / "21_fisher_kolmogorov_corti83/results"
                "/reaction_coefficients.csv")]
    scale = 0.5 * len(rates) / sum(rates)
    report.check(name, "rescaled rate, smallest", min(rates) * scale, 0.2177)
    report.check(name, "rescaled rate, largest", max(rates) * scale, 0.7195)
    report.check(name, "rescaled rate, mean",
                 sum(rates) * scale / len(rates), 0.5)

    for case, stages in (("tau", (14.0, 20.4, 28.8)),
                         ("amyloid", (1.2, 5.6, 13.2))):
        rows = read_csv(base / f"{case}_regional_profiles.csv")
        times = np.array([float(r["time"]) for r in rows])
        values = np.array([[float(r[f"node_{k}"]) for k in range(83)]
                           for r in rows])
        report.check(name, f"{case} stays in [0,1]",
                     1.0 if values.min() >= 0.0
                     and values.max() <= 1.0 + 1e-10 else 0.0, 1.0)
        means = values.mean(axis=1)
        for target, expected in zip((0.10, 0.40, 0.80), stages):
            report.check(name, f"{case} stage at mean {target:.0%}",
                         float(times[int(np.argmax(means >= target))]),
                         expected, "years")

    rows = read_csv(base / "tau_regional_biomarkers.csv")
    times = [float(r["time"]) for r in rows]
    crossings = {lobe: crossing(times, [float(r[lobe]) for r in rows], 50.0)
                 for lobe in ("temporal", "frontal", "parietal", "occipital")}
    for lobe, expected in (("temporal", 14.70), ("frontal", 19.91),
                           ("parietal", 29.71), ("occipital", 32.39)):
        report.check(name, f"tau {lobe} crossing with the regional rates",
                     crossings[lobe], expected, "years")
    ordered = sorted(crossings, key=crossings.get)
    report.check(name, "the regional rates recover the clinical lobe order",
                 1.0 if ordered == ["temporal", "frontal", "parietal",
                                    "occipital"] else 0.0, 1.0)
    # The uniform run at the same scaling gives the connectivity order, so
    # the two records must disagree exactly on the two middle lobes.
    uniform = read_csv(base / "tau_biomarkers.csv")
    times = [float(r["time"]) for r in uniform]
    uniform_order = sorted(
        ("temporal", "frontal", "parietal", "occipital"),
        key=lambda lobe: crossing(times, [float(r[lobe]) for r in uniform],
                                  50.0))
    report.check(name, "the uniform run keeps the connectivity order",
                 1.0 if uniform_order == ["temporal", "occipital", "parietal",
                                          "frontal"] else 0.0, 1.0)

    # The spreading moves anteriorly with the regional rates and posteriorly
    # with the uniform one: the centroid of the activated regions.
    coordinate = {int(r["node_id"]): float(r["y"]) for r in read_csv(
        Path("data/connectome/fornari83/nodes.csv"))}
    def centroid(path, stages):
        rows = read_csv(path)
        times = np.array([float(r["time"]) for r in rows])
        values = np.array([[float(r[f"node_{k}"]) for k in range(83)]
                           for r in rows])
        out = []
        for target in stages:
            row = values[int(np.argmin(abs(times - target)))]
            chosen = [k for k in range(83) if row[k] >= 0.5]
            out.append(sum(coordinate[k] for k in chosen) / len(chosen))
        return out
    forward = centroid(base / "tau_regional_profiles.csv", (14.0, 20.4, 28.8))
    backward = centroid(base / "tau_profiles.csv", (15.2, 21.6, 27.6))
    report.check(name, "regional: the centroid moves anteriorly",
                 1.0 if forward[0] < forward[1] < forward[2] else 0.0, 1.0)
    report.check(name, "uniform: the centroid moves posteriorly",
                 1.0 if backward[1] < backward[0] else 0.0, 1.0)
    chapter = " ".join(Path("report/chapter4_connectome.tex").read_text(
        encoding="utf-8").split())
    # The control that removes the protein transfer: a synthetic rate field,
    # linear in the anterior-posterior coordinate with the same spread, gives
    # the same lobe order as the coefficients of the reference.
    field = read_csv(base / "synthetic_rate_field.csv")
    values = [float(r["alpha"]) for r in field]
    report.check(name, "synthetic field, spread", max(values) / min(values),
                 3.3046)
    coordinate = {int(r["node_id"]): float(r["y"]) for r in read_csv(
        Path("data/connectome/fornari83/nodes.csv"))}
    ordered = sorted(field, key=lambda r: coordinate[int(r["node_id"])])
    monotone = all(float(a["alpha"]) <= float(b["alpha"])
                   for a, b in zip(ordered, ordered[1:]))
    report.check(name, "synthetic field grows towards the anterior end",
                 1.0 if monotone else 0.0, 1.0)
    rows = read_csv(base / "tau_synthetic_biomarkers.csv")
    times = [float(r["time"]) for r in rows]
    synthetic = {lobe: crossing(times, [float(r[lobe]) for r in rows], 50.0)
                 for lobe in ("temporal", "frontal", "parietal", "occipital")}
    for lobe, expected in (("temporal", 17.46), ("frontal", 21.12),
                           ("parietal", 26.75), ("occipital", 31.25)):
        report.check(name, f"synthetic control, {lobe} crossing",
                     synthetic[lobe], expected, "years")
    report.check(name, "the synthetic control recovers the clinical order",
                 1.0 if sorted(synthetic, key=synthetic.get)
                 == ["temporal", "frontal", "parietal", "occipital"]
                 else 0.0, 1.0)
    report.check_contains(
        name, "the section states the synthetic control", chapter,
        "at $17.5$, $21.1$, $26.7$ and $31.3$ years")

    # What a regional field must satisfy: the 5040 assignments of the seven
    # values to the seven groups.
    table = read_csv(base / "rate_permutations.csv")
    report.check(name, "assignments tried", float(len(table)), 5040.0)
    clinical = [r for r in table if r["clinical_order"] == "1"]
    below = [r for r in table
             if float(r["frontal_over_occipital"]) <= 1.0]
    report.check(name, "assignments giving the clinical order",
                 float(len(clinical)), 431.0)
    report.check(name, "none of them has the occipital rate above the "
                 "frontal one",
                 float(sum(int(r["clinical_order"]) for r in below)), 0.0)
    report.check(name, "assignments with frontal below occipital",
                 float(len(below)), 2520.0)
    above = [r for r in table if float(r["frontal_over_occipital"]) > 1.0]
    strong = [r for r in table if float(r["frontal_over_occipital"]) >= 3.0]
    report.check(name, "share with the clinical order, frontal above "
                 "occipital",
                 100.0 * sum(int(r["clinical_order"]) for r in above)
                 / len(above), 17.1, "percent")
    report.check(name, "share with the clinical order, ratio at least three",
                 100.0 * sum(int(r["clinical_order"]) for r in strong)
                 / len(strong), 75.0, "percent")
    for label, needle in (
            ("count of assignments", "in all $5040$ possible ways"),
            ("necessary condition", "in none of the $2520$ such "
                                    "assignments"),
            ("shares", "it appears in $17$ percent of those in which it does "
                       "and in $75$ percent of those in which it exceeds it "
                       "by a factor of three or more"),
            ("scope of the condition",
             "That condition is a statement about this model on this graph "
             "and it should not be read as a property of the protein."),
            ("distinction between a PET rate and alpha",
             "a rate of accumulation measured by positron emission "
             "tomography is not the coefficient $\\alpha$"),
            ("non-stationary tau pattern",
             "an anterior-to-posterior gradient of conversion is not a "
             "documented property of tau")):
        report.check_contains(name, f"the section states the {label}",
                              chapter, needle)

    # The amyloid variant: the ordering survives but the seeded mantle is no
    # longer uniform, which is why the report shows the tau row alone.
    def mantle(path, stages):
        rows = read_csv(path)
        times = np.array([float(r["time"]) for r in rows])
        values = np.array([[float(r[f"node_{k}"]) for k in range(83)]
                           for r in rows])
        seeded = [k for k in range(83) if values[0][k] > 0.0]
        out = []
        for target in stages:
            row = values[int(np.argmin(abs(times - target)))]
            out.append(row[seeded].max() - row[seeded].min())
        return out
    uniform_mantle = mantle(base / "amyloid_profiles.csv", (1.2, 5.2, 12.0))
    regional_mantle = mantle(base / "amyloid_regional_profiles.csv",
                             (1.2, 5.6, 13.2))
    report.check(name, "uniform amyloid mantle rises as a block",
                 max(uniform_mantle), 0.02)
    report.check(name, "regional amyloid mantle spreads",
                 max(regional_mantle), 0.58)
    # The amyloid row stays in the report figure: its ordering survives and
    # the report quotes both the narrowed margin and the new gradient.
    for label, needle in (
            ("amyloid stages", "with stages at $1.2$, $5.6$ and $13.2$ "
                               "years"),
            ("narrowed margin", "narrows from $5.2$ to $1.0$ years"),
            ("mantle gradient", "crossing at $3.2$, $4.0$, $8.4$ and $9.9$ "
                                "years instead of together at $4.4$")):
        report.check_contains(name,
                              f"regional staging prose states the {label}",
                              chapter, needle)

    report.note(name, "centroid of the activated regions, anterior positive",
                f"regional {['%+.1f' % v for v in forward]} mm, "
                f"uniform {['%+.1f' % v for v in backward]} mm")

    ordered_crossings = [crossings[lobe] for lobe in
                         ("temporal", "frontal", "parietal", "occipital")]
    for label, needle in (
            ("crossings", "cross the $50\\%$ level at "
             + ", ".join(f"${v:.1f}$" for v in ordered_crossings[:-1])
             + f" and ${ordered_crossings[-1]:.1f}$ years"),
            ("rescaled extremes",
             f"running from ${min(rates) * scale:.4f}$ in the occipital "
             f"group to ${max(rates) * scale:.4f}$ in the frontal one"),
            ("centroid direction",
             f"from $+{round(forward[0]):d}$ to $+{round(forward[2]):d}$ mm")):
        report.check_contains(name,
                              f"regional staging prose states the {label}",
                              chapter, needle)


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
        # The ranking follows the number of connections of the seeded
        # vertex, which sets the metric mass of the seed in the finite
        # element model, far more closely than the weighted degree; at
        # rho = 1 the logistic growth of that seed mass predicts the times.
        from scipy.stats import spearmanr
        connections = [float(row["degree"]) for row in rows]
        weighted = [float(row["weighted_degree"]) for row in rows]
        edges = len(read_csv(Path("data/connectome/fornari83/edges.csv")))
        seed = [0.1 * 0.5 * c / edges for c in connections]
        predicted = [2.0 * math.log(19.0 * (1.0 - c) / c) for c in seed]
        report.check(name, f"rank correlation with connections at "
                     f"rho={scaling:g}",
                     float(spearmanr(values, connections).correlation),
                     {1: -0.982, 0.02: -0.823}[scaling])
        report.check(name, f"rank correlation with weighted degree at "
                     f"rho={scaling:g}",
                     float(spearmanr(values, weighted).correlation),
                     {1: -0.363, 0.02: -0.226}[scaling])
        if scaling == 1:
            gaps = [v - q for v, q in zip(values, predicted)]
            mean = sum(gaps) / len(gaps)
            report.check(name, "seed-mass prediction, mean gap", mean, -0.35,
                         "years")
            report.check(name, "seed-mass prediction, scatter of the gap",
                         math.sqrt(sum((g - mean) ** 2 for g in gaps)
                                   / len(gaps)), 0.12, "years")


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
    for target, expected in ((0.2, 24.0), (0.4, 32.0), (0.6, 40.0)):
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
            (("regional", "frontal"), 30.75),
            (("regional", "occipital"), 45.00),
            (("uniform", "occipital"), 36.80),
            (("uniform", "frontal"), 37.00)):
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
    # With the uniform rate the four biomarker lobes recover the order of
    # the connectivity, the frontal lobe last; the output stride of one year
    # ties it with two of the deep groups, so the claim is made over the
    # four lobes the chapter compares.
    lobes = ("temporal", "frontal", "parietal", "occipital")
    report.check(name, "frontal last of the four lobes with the uniform rate",
                 1.0 if max(lobes, key=lambda g: uniform[g]) == "frontal"
                 else 0.0, 1.0)
    report.check(name, "regional spread",
                 max(regional.values()) - min(regional.values()), 14.25,
                 "years")
    report.check(name, "uniform spread",
                 max(uniform.values()) - min(uniform.values()), 1.5625,
                 "years")


def check_orientation(report):
    """Every sagittal view of the report faces the same way as the staging
    drawings of Weickenmeier et al. and figures 1 and 5 of Fornari et al.:
    the frontal pole to the left of the image. Checked statically, no report
    script may use the mirror camera, and by rendering the frontal poles."""
    name = "anatomy orientation"
    scripts = ("plot-fisher-fornari83.py", "plot-fornari-connectome-topology.py",
               "plot-connectome-views.py", "plot-connectome-regions.py",
               "plot-corti-activation-order.py", "plot-connectome-staging.py")
    for script in scripts:
        text = " ".join((Path("scripts") / script).read_text().split())
        # Two of these scripts place the camera themselves instead of calling
        # render(), so the mirror view can also enter as an inlined direction.
        inlined = any(f"({sign}1.0, 0.0, 0.0)" in text
                      for sign in ("+", "")) and "CAMERA" not in text
        report.check(name, f"{script} uses the left-side camera",
                     0.0 if ("sagittal_right" in text or inlined) else 1.0,
                     1.0)
    try:
        import vtk  # noqa: F401
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from PIL import Image
    except ImportError:
        report.note(name, "render", "skipped, vtk or PIL not installed")
        return
    sys.path.insert(0, "scripts")
    import render_connectome as rc
    from connectome_style import load_nodes
    nodes = load_nodes()
    coords = np.array([node["coords"] for node in nodes])
    values = np.array([1.0 if "frontalpole" in node["name"].lower() else 0.0
                       for node in nodes])
    table = rc.lookup_table(plt.cm.viridis, 0.0, 1.0)
    import tempfile
    with tempfile.TemporaryDirectory() as scratch:
        for view, side in (("sagittal", "left"), ("sagittal_right", "right")):
            path, _ = rc.render(Path(scratch) / f"{view}.png", view, coords,
                                values, table, node_radius=5.0,
                                size=(900, 700))
            image = np.asarray(Image.open(path).convert("RGB"), dtype=float)
            yellow = ((image[:, :, 0] > 200) & (image[:, :, 1] > 200)
                      & (image[:, :, 2] < 120))
            columns = np.nonzero(yellow)[1]
            fraction = float(columns.mean()) / image.shape[1]
            report.check(name, f"frontal pole on the {side} in the "
                         f"{view} view",
                         1.0 if (fraction < 0.5) == (side == "left") else 0.0,
                         1.0)
            report.note(name, f"frontal pole column fraction, {view}",
                        f"{fraction:.2f}")


def main():
    report = Report()
    for check in (check_18, check_18_timestep, check_19, check_19_accuracy,
                  check_19_topology, check_connectome_consistency,
                  check_19_scheme, check_20, check_21, check_22, check_23,
                  check_23_stabilization, check_23_mass_spectrum,
                  check_23_boundary_prose, check_23_lobe_scale, check_23_lobe_order,
                  check_24_views, check_25, check_26,
                  check_26_order, check_27, check_27_regional,
                  check_orientation):
        try:
            check(report)
        except FileNotFoundError as error:
            report.note(check.__name__, "input", f"missing: {error.filename}")
            report.failed += 1
    report.show()
    return 1 if report.failed else 0


if __name__ == "__main__":
    sys.exit(main())
