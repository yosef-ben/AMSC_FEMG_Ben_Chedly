#!/usr/bin/env python3

"""Is the connectivity-driven lobe order robust to the coarse-graining?

The chapter derives the uniform-rate activation order of the four lobes,
temporal, occipital, parietal, frontal, on the 83-region graph obtained by
filtering the fine Budapest graph at five occurrences and aggregating the
parcels by FreeSurfer parent region, with parallel connections summed. This
study repeats the reference tau run on the FINE graph, the 1015 parcels with
the same five-occurrence filter and no aggregation, and asks whether the same
order comes out. Everything else is held fixed and nothing is tuned: the same
model (metric-graph FEM, lumped mass, one P1 element per unit-length
connection, backward Euler with Newton), the same alpha = 0.5, the same
rho = 0.005, the same dt = 0.4 over 80 years, the seed c0 = 0.1 on the
parcels whose parent region is an entorhinal cortex, and the lobes assigned
to the parcels through the parent-region names with the classify rule of the
report. The lobe biomarkers are aggregated only after the simulation, in two
ways: the mean over the parcel curves of the lobe, and the mean over the
regions of the lobe of each region's parcel mean, which matches the
equal-weight regional biomarker of the report. The integrator is validated
first by reproducing the stored 83-region run to solver tolerance.

Internal robustness check only: it does not touch the report.
"""

import csv
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from lobe_scale import classify  # noqa: E402  (the report's lobe rule)

GRAPHML = ROOT / "data/connectome/source/budapest_all_20k.graphml"
STORED = ROOT / "benchmarks/27_connectome_seeding_patterns/results"
RESULTS = Path(__file__).resolve().parent / "results"
NS = {"g": "http://graphml.graphdrawing.org/xmlns"}
MIN_OCCURRENCES = 5
ALPHA = 0.5
RHO = 0.005
DT = 0.4
FINAL_TIME = 80.0
LOBES = ("temporal", "occipital", "parietal", "frontal")


def parent_name(node):
    name = re.sub(r"_[0-9]+$", "", node["dn_name"])
    if name.startswith("rh."):
        return "ctx-rh-" + name[3:]
    if name.startswith("lh."):
        return "ctx-lh-" + name[3:]
    return name


def read_fine_graph():
    root = ET.parse(GRAPHML).getroot()
    keys = {item.attrib["id"]: item.attrib.get("attr.name", item.attrib["id"])
            for item in root.findall("g:key", NS)}
    graph = root.find("g:graph", NS)
    nodes = {}
    for node in graph.findall("g:node", NS):
        nodes[node.attrib["id"]] = {
            keys[data.attrib["key"]]: data.text or ""
            for data in node.findall("g:data", NS)}
    edges = []
    for edge in graph.findall("g:edge", NS):
        values = {keys[data.attrib["key"]]: data.text or ""
                  for data in edge.findall("g:data", NS)}
        if int(values["occurences"]) < MIN_OCCURRENCES:
            continue
        weight = (float(values["fiber_count_median"])
                  / float(values["fiber_length_median"]))
        edges.append((edge.attrib["source"], edge.attrib["target"], weight))
    return nodes, edges


def backward_euler_lumped(laplacian, lumped, alpha, initial, dt, final_time):
    """The be_lumped scheme of the library: M_L (c - c_old) + dt rho L c
    = dt M_L alpha c (1 - c), solved by Newton at every step."""
    steps = int(round(final_time / dt))
    state = initial.copy()
    history = [state.copy()]
    size = state.size
    for _ in range(steps):
        old = state.copy()
        iterate = state.copy()
        scale = max(1.0, np.linalg.norm(lumped * old))
        for _ in range(60):
            residual = (lumped * (iterate - old) + dt * (laplacian @ iterate)
                        - dt * lumped * alpha * iterate * (1.0 - iterate))
            if np.linalg.norm(residual) <= 1.0e-11 * scale:
                break
            diagonal = lumped * (1.0 - dt * alpha * (1.0 - 2.0 * iterate))
            jacobian = csr_matrix(
                (dt * laplacian).tocsr()
                + csr_matrix((diagonal, (np.arange(size), np.arange(size)))))
            iterate = iterate + splu(jacobian.tocsc()).solve(-residual)
        state = iterate
        history.append(state.copy())
    times = np.arange(steps + 1) * dt
    return times, np.array(history)


def crossing(times, series, level=0.5):
    for k in range(1, len(series)):
        if series[k - 1] < level <= series[k]:
            return times[k - 1] + (level - series[k - 1]) / (
                series[k] - series[k - 1]) * (times[k] - times[k - 1])
    return float("nan")


def build_operator(n, edges):
    laplacian = np.zeros((n, n))
    lumped = np.zeros(n)
    for i, j, w in edges:
        laplacian[i, i] += w
        laplacian[j, j] += w
        laplacian[i, j] -= w
        laplacian[j, i] -= w
        lumped[i] += 0.5
        lumped[j] += 0.5
    return csr_matrix(RHO * laplacian), lumped


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)

    # --- validation: reproduce the stored 83-region run -------------------
    names = {int(r["node_id"]): r["name"] for r in csv.DictReader(
        open(ROOT / "data/connectome/fornari83/nodes.csv"))}
    coarse_edges = [(int(e["source"]), int(e["target"]),
                     float(e["connectivity_weight"]))
                    for e in csv.DictReader(
                        open(ROOT / "data/connectome/fornari83/edges.csv"))]
    laplacian, lumped = build_operator(83, coarse_edges)
    initial = np.zeros(83)
    for k, name in names.items():
        if "entorhinal" in name.lower():
            initial[k] = 0.1
    times, history = backward_euler_lumped(
        laplacian, lumped, ALPHA, initial, DT, FINAL_TIME)
    lobe_of = {k: classify(name) for k, name in names.items()}
    stored = list(csv.DictReader(open(STORED / "tau_biomarkers.csv")))
    stored_times = [float(r["time"]) for r in stored]
    print("validation against the stored 83-region run:")
    coarse_crossings = {}
    for lobe in LOBES:
        members = [k for k in names if lobe_of[k] == lobe]
        ours = crossing(times, history[:, members].mean(axis=1))
        reference = crossing(stored_times,
                             [float(r[lobe]) / 100.0 for r in stored])
        coarse_crossings[lobe] = reference
        print(f"  {lobe:9s} integrator {ours:7.3f}  stored {reference:7.3f}"
              f"  difference {abs(ours - reference):.2e}")
        assert abs(ours - reference) < 1.0e-6, "integrator validation failed"

    # --- the fine graph ---------------------------------------------------
    nodes, raw_edges = read_fine_graph()
    index = {key: k for k, key in enumerate(sorted(nodes))}
    region = {index[key]: parent_name(nodes[key]) for key in index}
    lobe = {k: classify(region[k]) for k in region}
    fine_edges = [(index[a], index[b], w) for a, b, w in raw_edges]
    n = len(index)
    print(f"\nfine graph: {n} parcels, {len(fine_edges)} connections "
          f"after the {MIN_OCCURRENCES}-occurrence filter")
    laplacian, lumped = build_operator(n, fine_edges)
    initial = np.zeros(n)
    seeds = [k for k in region if "entorhinal" in region[k].lower()]
    initial[seeds] = 0.1
    print(f"seed: {len(seeds)} parcels of the entorhinal regions")
    times, history = backward_euler_lumped(
        laplacian, lumped, ALPHA, initial, DT, FINAL_TIME)

    # --- post-hoc aggregation --------------------------------------------
    with open(RESULTS / "fine_lobe_crossings.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["lobe", "parcel_mean_crossing",
                         "region_mean_crossing", "coarse_crossing"])
        rows = {}
        for lobe_name in LOBES:
            parcels = [k for k in region if lobe[k] == lobe_name]
            parcel_mean = crossing(times, history[:, parcels].mean(axis=1))
            regions = sorted({region[k] for k in parcels})
            per_region = np.column_stack([
                history[:, [k for k in parcels if region[k] == r]].mean(axis=1)
                for r in regions])
            region_mean = crossing(times, per_region.mean(axis=1))
            rows[lobe_name] = (parcel_mean, region_mean)
            writer.writerow([lobe_name, f"{parcel_mean:.4f}",
                             f"{region_mean:.4f}",
                             f"{coarse_crossings[lobe_name]:.4f}"])

    print("\nlobe crossings of the mean curves [years]:")
    print(f"  {'lobe':9s} {'fine, parcel mean':>18s} "
          f"{'fine, region mean':>18s} {'coarse (stored)':>16s}")
    for lobe_name in LOBES:
        parcel_mean, region_mean = rows[lobe_name]
        print(f"  {lobe_name:9s} {parcel_mean:18.2f} {region_mean:18.2f}"
              f" {coarse_crossings[lobe_name]:16.2f}")
    for label, column in (("parcel-mean", 0), ("region-mean", 1)):
        order = sorted(LOBES, key=lambda l: rows[l][column])
        values = [rows[l][column] for l in order]
        print(f"  fine order ({label}): {', '.join(order)}; "
              f"separation {max(values) - min(values):.2f} years")
    coarse_order = sorted(LOBES, key=coarse_crossings.get)
    print(f"  coarse order: {', '.join(coarse_order)}; separation "
          f"{max(coarse_crossings.values()) - min(coarse_crossings.values()):.2f}"
          " years")

    # --- per-region correlation between the two levels --------------------
    coarse_profiles = list(csv.DictReader(open(STORED / "tau_profiles.csv")))
    coarse_t = [float(r["time"]) for r in coarse_profiles]
    fine_region_time = {}
    coarse_region_time = {}
    for k, name in names.items():
        parcels = [p for p in region if region[p] == name]
        if not parcels:
            continue
        fine_region_time[name] = crossing(
            times, history[:, parcels].mean(axis=1))
        coarse_region_time[name] = crossing(
            coarse_t, [float(r[f"node_{k}"]) for r in coarse_profiles])
    shared = [name for name in fine_region_time
              if not np.isnan(fine_region_time[name])
              and not np.isnan(coarse_region_time[name])]
    rho = spearmanr([fine_region_time[name] for name in shared],
                    [coarse_region_time[name] for name in shared]).correlation
    print(f"\nper-region activation times, fine against coarse: "
          f"Spearman {rho:.3f} over {len(shared)} regions")
    with open(RESULTS / "fine_region_times.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["region", "fine_years", "coarse_years"])
        for name in sorted(shared):
            writer.writerow([name, f"{fine_region_time[name]:.4f}",
                             f"{coarse_region_time[name]:.4f}"])


if __name__ == "__main__":
    main()
