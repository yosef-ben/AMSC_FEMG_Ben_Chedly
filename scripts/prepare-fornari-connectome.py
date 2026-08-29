#!/usr/bin/env python3

"""Aggregate the Budapest v3 connectome to the 83-region FreeSurfer graph.

The weight of a fine connection is its median fibre count divided by its
median fibre length, and parallel connections between the same two regions
add, as conductances in parallel. Although Fornari et al. describe the fibre
quantities as cohort means, it is the median fields of the public file that
reproduce every published graph statistic to printed precision, including
which regions attain the extremes; the mean fields do not. The per-edge
``electrical_connectivity_median`` field, the median of the per-subject
ratios, is close but reproduces none of the published values exactly.

The five-subject occurrence threshold is not stated in the reference either;
it is the only integer that reproduces the published fine graph of 1015
vertices and 37477 edges (4 keeps 40895 edges, 6 keeps 34718).
"""

import argparse
import csv
import json
import math
import re
import statistics
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


NS = {"g": "http://graphml.graphdrawing.org/xmlns"}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--graphml",
        type=Path,
        default=Path("data/connectome/source/budapest_all_20k.graphml"),
    )
    parser.add_argument(
        "--anatomy",
        type=Path,
        default=Path("data/connectome/anatomy/nodes.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/connectome/fornari83"),
    )
    parser.add_argument(
        "--minimum-occurrences",
        type=int,
        default=5,
        help="Keep fibres found in at least this many subjects (default: 5).",
    )
    return parser.parse_args()


def read_graphml(path):
    root = ET.parse(path).getroot()
    keys = {
        item.attrib["id"]: item.attrib.get("attr.name", item.attrib["id"])
        for item in root.findall("g:key", NS)
    }
    graph = root.find("g:graph", NS)
    if graph is None or graph.attrib.get("edgedefault") != "undirected":
        raise ValueError("Expected one undirected GraphML graph")

    nodes = {}
    for node in graph.findall("g:node", NS):
        nodes[node.attrib["id"]] = {
            keys[data.attrib["key"]]: data.text or ""
            for data in node.findall("g:data", NS)
        }

    edges = []
    for edge in graph.findall("g:edge", NS):
        values = {
            keys[data.attrib["key"]]: data.text or ""
            for data in edge.findall("g:data", NS)
        }
        edges.append((edge.attrib["source"], edge.attrib["target"], values))
    return nodes, edges


def parent_name(node):
    name = re.sub(r"_[0-9]+$", "", node["dn_name"])
    if name.startswith("rh."):
        return "ctx-rh-" + name[3:]
    if name.startswith("lh."):
        return "ctx-lh-" + name[3:]
    return name


def parent_key(node):
    name = parent_name(node)
    hemisphere = "midline" if name == "Brain-Stem" else node["dn_hemisphere"]
    return hemisphere, name


def read_anatomy(path):
    raw = json.loads(path.read_text(encoding="utf-8"))
    regions = []
    for source_id, values in sorted(raw.items(), key=lambda item: int(item[0])):
        name = values["dn_fsname"]
        hemisphere = "midline" if name == "Brain-Stem" else values["dn_hemisphere"]
        regions.append(
            {
                "source_id": int(source_id),
                "key": (hemisphere, name),
                "name": name,
                "hemisphere": hemisphere,
                "region": values["dn_region"],
                "x": float(values["pial_x"]),
                "y": float(values["pial_y"]),
                "z": float(values["pial_z"]),
            }
        )
    return regions


def summary(values):
    return {
        "minimum": min(values),
        "maximum": max(values),
        "mean": statistics.mean(values),
    }


def aggregate(nodes, edges, regions, minimum_occurrences):
    known_regions = {region["key"] for region in regions}
    fine_edges = []
    grouped = defaultdict(list)
    intra_region_edges = 0

    for source, target, values in edges:
        if int(values["occurences"]) < minimum_occurrences:
            continue
        fine_edges.append((source, target, values))
        source_region = parent_key(nodes[source])
        target_region = parent_key(nodes[target])
        if source_region not in known_regions or target_region not in known_regions:
            raise ValueError(f"Fine node cannot be mapped to atlas: {source_region} {target_region}")
        if source_region == target_region:
            intra_region_edges += 1
            continue
        pair = tuple(sorted((source_region, target_region)))
        grouped[pair].append(values)

    coarse_edges = []
    for pair, values in sorted(grouped.items()):
        fibre_number = sum(float(value["fiber_count_median"]) for value in values)
        fibre_length = statistics.mean(
            float(value["fiber_length_median"]) for value in values
        )
        connectivity = sum(
            float(value["fiber_count_median"]) / float(value["fiber_length_median"])
            for value in values
        )
        coarse_edges.append(
            {
                "source_key": pair[0],
                "target_key": pair[1],
                "fine_edges": len(values),
                "fibre_number": fibre_number,
                "fibre_length": fibre_length,
                "connectivity": connectivity,
            }
        )

    return fine_edges, intra_region_edges, coarse_edges


def write_outputs(output_dir, regions, coarse_edges, metadata):
    output_dir.mkdir(parents=True, exist_ok=True)
    region_id = {region["key"]: index for index, region in enumerate(regions)}

    with (output_dir / "nodes.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "node_id", "source_id", "name", "hemisphere", "region", "x", "y", "z"
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for index, region in enumerate(regions):
            writer.writerow(
                {
                    "node_id": index,
                    "source_id": region["source_id"],
                    "name": region["name"],
                    "hemisphere": region["hemisphere"],
                    "region": region["region"],
                    "x": f"{region['x']:.12g}",
                    "y": f"{region['y']:.12g}",
                    "z": f"{region['z']:.12g}",
                }
            )

    with (output_dir / "edges.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "edge_id", "source", "target", "fine_edges", "fibre_number",
            "fibre_length_mm", "connectivity_weight"
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for index, edge in enumerate(coarse_edges):
            writer.writerow(
                {
                    "edge_id": index,
                    "source": region_id[edge["source_key"]],
                    "target": region_id[edge["target_key"]],
                    "fine_edges": edge["fine_edges"],
                    "fibre_number": f"{edge['fibre_number']:.12g}",
                    "fibre_length_mm": f"{edge['fibre_length']:.12g}",
                    "connectivity_weight": f"{edge['connectivity']:.12g}",
                }
            )

    for cells_per_edge in (1, 2, 4, 8, 16, 32, 64):
        graph_path = output_dir / f"graph_fem_{cells_per_edge}.txt"
        with graph_path.open("w", encoding="ascii") as stream:
            stream.write(f"{len(regions)} {len(coarse_edges)}\n")
            for edge in coarse_edges:
                source = region_id[edge["source_key"]]
                target = region_id[edge["target_key"]]
                stream.write(f"{source} {target} 1 {cells_per_edge}\n")
    (output_dir / "graph_fem.txt").write_text(
        (output_dir / "graph_fem_1.txt").read_text(encoding="ascii"),
        encoding="ascii",
    )

    (output_dir / "summary.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def main():
    args = arguments()
    nodes, edges = read_graphml(args.graphml)
    regions = read_anatomy(args.anatomy)
    fine_edges, intra_region_edges, coarse_edges = aggregate(
        nodes, edges, regions, args.minimum_occurrences
    )

    fibre_numbers = [edge["fibre_number"] for edge in coarse_edges]
    fibre_lengths = [edge["fibre_length"] for edge in coarse_edges]
    connectivities = [edge["connectivity"] for edge in coarse_edges]
    degrees = [0] * len(regions)
    weighted_degrees = [0.0] * len(regions)
    region_id = {region["key"]: index for index, region in enumerate(regions)}
    for edge in coarse_edges:
        source = region_id[edge["source_key"]]
        target = region_id[edge["target_key"]]
        degrees[source] += 1
        degrees[target] += 1
        weighted_degrees[source] += edge["connectivity"]
        weighted_degrees[target] += edge["connectivity"]

    # Statistics of the retained fine connections themselves. The reference
    # prints a fibre-length range of 11.3-136.8 mm: the lower bound and the
    # mean are reproduced by the region-to-region connections above, the
    # upper bound by the longest retained fine connection recorded here.
    fine_lengths = [float(values["fiber_length_median"])
                    for _, _, values in fine_edges]
    fine_counts = [float(values["fiber_count_median"])
                   for _, _, values in fine_edges]

    metadata = {
        "source": str(args.graphml),
        "minimum_occurrences": args.minimum_occurrences,
        "fine_graph": {
            "vertices": len(nodes),
            "retained_edges": len(fine_edges),
            "fibre_length_mm": summary(fine_lengths),
            "fibre_number": summary(fine_counts),
        },
        "aggregation": {
            "atlas": "FreeSurfer 83-region parcellation",
            "vertices": len(regions),
            "intra_region_edges_removed": intra_region_edges,
            "edges": len(coarse_edges),
            "degree": summary(degrees),
            "weighted_degree": summary(weighted_degrees),
            "fibre_number": summary(fibre_numbers),
            "fibre_length_mm": summary(fibre_lengths),
            "connectivity_weight": summary(connectivities),
        },
        "paper_reference": {
            "fine_vertices": 1015,
            "fine_edges": 37477,
            "coarse_vertices": 83,
            "coarse_edges": 1130,
            "degree": {"minimum": 6, "maximum": 48},
            "fibre_number": 40.2,
            "fibre_number_range": {"minimum": 1, "maximum": 596},
            "fibre_length_mm": 38.4,
            "fibre_length_range_mm": {"minimum": 11.3, "maximum": 136.8},
            "connectivity_weight": {
                "minimum": 0.01, "maximum": 35.32, "mean": 1.57
            },
            "weighted_degree": {
                "minimum": 2.1, "maximum": 127.6, "mean": 42.8
            },
        },
    }

    expected = (1015, 37477, 83, 1130)
    actual = (len(nodes), len(fine_edges), len(regions), len(coarse_edges))
    if actual != expected:
        raise ValueError(f"Unexpected topology {actual}; expected {expected}")
    if not math.isclose(statistics.mean(fibre_numbers), 40.2, abs_tol=0.1):
        raise ValueError("Aggregated fibre-number statistics do not match the paper")
    if not math.isclose(statistics.mean(fibre_lengths), 38.4, abs_tol=0.1):
        raise ValueError("Aggregated fibre-length statistics do not match the paper")

    write_outputs(args.output_dir, regions, coarse_edges, metadata)
    print(f"Fine graph:   {len(nodes)} vertices, {len(fine_edges)} edges")
    print(f"Coarse graph: {len(regions)} vertices, {len(coarse_edges)} edges")
    print(f"Mean fibres:  {statistics.mean(fibre_numbers):.6f}")
    print(f"Mean length:  {statistics.mean(fibre_lengths):.6f} mm")


if __name__ == "__main__":
    main()
