#!/usr/bin/env python3

"""Compare published and reproducible 83-region connectivity statistics."""

import argparse
import csv
import importlib.util
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--median-edges", type=Path,
        default=Path("data/connectome/fornari83/edges.csv"))
    parser.add_argument("--mean-graphml", type=Path,
        default=Path("data/connectome/source/budapest_all_20k_mean_electrical.graphml"))
    parser.add_argument("--anatomy", type=Path,
        default=Path("data/connectome/anatomy/nodes.json"))
    parser.add_argument("--output", type=Path,
        default=Path("benchmarks/19_fisher_kolmogorov_fornari83/results/connectome_weight_audit.csv"))
    return parser.parse_args()


def load_preprocessor():
    path = Path(__file__).with_name("prepare-fornari-connectome.py")
    spec = importlib.util.spec_from_file_location("fornari_preprocessor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def statistics_row(name, node_count, edge_rows):
    weights = [weight for _, _, weight in edge_rows]
    degrees = [0.0] * node_count
    laplacian = np.zeros((node_count, node_count))
    for source, target, weight in edge_rows:
        degrees[source] += weight
        degrees[target] += weight
        laplacian[source, source] += weight
        laplacian[target, target] += weight
        laplacian[source, target] -= weight
        laplacian[target, source] -= weight
    eigenvalues = np.linalg.eigvalsh(laplacian)
    return {
        "variant": name, "nodes": node_count, "edges": len(edge_rows),
        "adjacency_min": min(weights), "adjacency_max": max(weights),
        "adjacency_mean": statistics.mean(weights),
        "weighted_degree_min": min(degrees),
        "weighted_degree_max": max(degrees),
        "weighted_degree_mean": statistics.mean(degrees),
        "lambda_1": eigenvalues[1],
    }


def read_median_edges(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return [(int(row["source"]), int(row["target"]),
             float(row["connectivity_weight"])) for row in rows]


def aggregate_mean_edges(graphml, anatomy):
    preprocessor = load_preprocessor()
    nodes, edges = preprocessor.read_graphml(graphml)
    regions = preprocessor.read_anatomy(anatomy)
    region_id = {region["key"]: index for index, region in enumerate(regions)}
    grouped = defaultdict(float)
    for source, target, values in edges:
        source_region = preprocessor.parent_key(nodes[source])
        target_region = preprocessor.parent_key(nodes[target])
        if source_region == target_region:
            continue
        pair = tuple(sorted((source_region, target_region)))
        grouped[pair] += float(values["number_of_fiber_per_fiber_length_mean"])
    coarse_edges = [(region_id[pair[0]], region_id[pair[1]], weight)
                    for pair, weight in sorted(grouped.items())]
    return len(regions), coarse_edges


def main():
    args = arguments()
    median_edges = read_median_edges(args.median_edges)
    node_count, mean_edges = aggregate_mean_edges(args.mean_graphml, args.anatomy)
    rows = [
        {
            "variant": "Fornari et al. (reported)", "nodes": 83, "edges": 1130,
            "adjacency_min": 0.01, "adjacency_max": 35.32,
            "adjacency_mean": 1.57, "weighted_degree_min": 2.1,
            "weighted_degree_max": 127.6, "weighted_degree_mean": 42.8,
            "lambda_1": "",
        },
        statistics_row("public median aggregation", node_count, median_edges),
        statistics_row("official mean-mode aggregation", node_count, mean_edges),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['variant']}: adjacency mean={row['adjacency_mean']}, "
              f"degree mean={row['weighted_degree_mean']}, "
              f"lambda_1={row['lambda_1']}")


if __name__ == "__main__":
    main()
