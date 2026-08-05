#!/usr/bin/env python3

"""Validate and convert the Budapest Reference Connectome for FEMG."""

import argparse
import csv
import json
import math
import statistics
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


NS = {"g": "http://graphml.graphdrawing.org/xmlns"}
GRAPHML_WEIGHT = "number_of_fiber_per_fiber_length_mean"
CSV_WEIGHT = "edge weight(med nof)"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graphml", type=Path,
                        default=Path("data/budapest_connectome_3.0_209_0_median.graphml"))
    parser.add_argument("--csv", type=Path,
                        default=Path("data/budapest_connectome_3.0_209_0_median.csv"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("data/connectome"))
    parser.add_argument("--edge-length", type=float, default=1.0)
    parser.add_argument("--cells-per-edge", type=int, default=1)
    return parser.parse_args()


def id_key(value):
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)


def edge_key(source, target):
    return tuple(sorted((source, target), key=id_key))


def read_graphml(path):
    root = ET.parse(path).getroot()
    names = {item.attrib["id"]: item.attrib.get("attr.name", item.attrib["id"])
             for item in root.findall("g:key", NS)}
    graph = root.find("g:graph", NS)
    if graph is None or graph.attrib.get("edgedefault") != "undirected":
        raise ValueError("Expected one undirected GraphML graph")

    nodes = {}
    for node in graph.findall("g:node", NS):
        nodes[node.attrib["id"]] = {
            names[data.attrib["key"]]: data.text or ""
            for data in node.findall("g:data", NS)
        }

    edges = []
    for edge in graph.findall("g:edge", NS):
        values = {names[data.attrib["key"]]: data.text or ""
                  for data in edge.findall("g:data", NS)}
        source, target = edge.attrib["source"], edge.attrib["target"]
        if source not in nodes or target not in nodes:
            raise ValueError("GraphML edge references an unknown node")
        edges.append({"source_original": source, "target_original": target,
                      "weight": float(values[GRAPHML_WEIGHT])})
    return nodes, edges


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter=";")
        expected = {"id node1", "id node2", "parent id node1", "parent id node2",
                    "parent name node1", "parent name node2", CSV_WEIGHT}
        if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
            raise ValueError("Unexpected Budapest CSV header")
        return list(reader)


def validate_and_merge(nodes, graphml_edges, csv_rows):
    graphml_values = {
        edge_key(edge["source_original"], edge["target_original"]): edge["weight"]
        for edge in graphml_edges
    }
    csv_values = {}
    for row in csv_rows:
        source, target = row["id node1"], row["id node2"]
        key = edge_key(source, target)
        if key in csv_values:
            raise ValueError(f"Duplicate CSV edge {key}")
        csv_values[key] = float(row[CSV_WEIGHT])
        for suffix, node_id in (("1", source), ("2", target)):
            if node_id not in nodes:
                raise ValueError(f"CSV node {node_id} is missing from GraphML")
            nodes[node_id].update({
                "parent_id": row[f"parent id node{suffix}"],
                "parent_name": row[f"parent name node{suffix}"],
            })

    if len(graphml_values) != len(graphml_edges):
        raise ValueError("Duplicate undirected edge in GraphML")
    if graphml_values.keys() != csv_values.keys():
        raise ValueError("GraphML and CSV contain different edge sets")
    if any(not math.isclose(graphml_values[key], csv_values[key], abs_tol=1e-12)
           for key in graphml_values):
        raise ValueError("GraphML and CSV contain different edge weights")


def connected_components(nodes, edges):
    adjacency = {node_id: set() for node_id in nodes}
    for edge in edges:
        source, target = edge["source_original"], edge["target_original"]
        adjacency[source].add(target)
        adjacency[target].add(source)
    output, visited = [], set()
    for start in sorted(nodes, key=id_key):
        if start in visited:
            continue
        stack, current = [start], []
        visited.add(start)
        while stack:
            node = stack.pop()
            current.append(node)
            for neighbor in adjacency[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        output.append(current)
    output.sort(key=lambda item: -len(item))
    return output, adjacency


def percentile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def weight_summary(values):
    return {"minimum": min(values), "first_quartile": percentile(values, 0.25),
            "median": statistics.median(values),
            "third_quartile": percentile(values, 0.75), "maximum": max(values),
            "mean": statistics.mean(values)}


def make_summary(nodes, edges, components, adjacency, selected, retained):
    degree = Counter()
    for edge in retained:
        degree[edge["source_original"]] += 1
        degree[edge["target_original"]] += 1
    degrees = [degree[node] for node in selected]
    return {
        "source": {
            "vertices": len(nodes), "edges": len(edges),
            "connected_components": len(components),
            "component_sizes": [len(item) for item in components],
            "isolated_vertices": sum(not adjacency[node] for node in nodes),
            "self_loops": sum(edge["source_original"] == edge["target_original"]
                              for edge in edges),
            "regions": dict(sorted(Counter(x.get("dn_region", "")
                                            for x in nodes.values()).items())),
            "hemispheres": dict(sorted(Counter(x.get("dn_hemisphere", "")
                                                for x in nodes.values()).items())),
            "weight": weight_summary([edge["weight"] for edge in edges]),
        },
        "processed": {
            "selection": "largest connected component", "vertices": len(selected),
            "edges_before_self_loop_removal": sum(
                edge["source_original"] in selected and edge["target_original"] in selected
                for edge in edges),
            "removed_self_loops": sum(
                edge["source_original"] == edge["target_original"]
                and edge["source_original"] in selected for edge in edges),
            "edges": len(retained),
            "degree": {"minimum": min(degrees), "maximum": max(degrees),
                       "mean": statistics.mean(degrees)},
            "weight": weight_summary([edge["weight"] for edge in retained]),
        },
        "validation": {"graphml_csv_topology_match": True,
                       "graphml_csv_weights_match": True},
    }


def write_outputs(output_dir, nodes, selected, retained, summary, length, cells):
    output_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(selected, key=id_key)
    local_id = {node: index for index, node in enumerate(ordered)}
    coordinates = {
        node: (math.cos(2.0 * math.pi * index / len(ordered)),
               math.sin(2.0 * math.pi * index / len(ordered)))
        for index, node in enumerate(ordered)
    }

    with (output_dir / "budapest_lcc_nodes.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        fields = ["local_id", "original_id", "name", "freesurfer_name", "region",
                  "hemisphere", "parent_id", "parent_name", "layout_x", "layout_y"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for node in ordered:
            data, (x, y) = nodes[node], coordinates[node]
            writer.writerow({
                "local_id": local_id[node], "original_id": node,
                "name": data.get("dn_name", ""),
                "freesurfer_name": data.get("dn_fsname", ""),
                "region": data.get("dn_region", ""),
                "hemisphere": data.get("dn_hemisphere", ""),
                "parent_id": data.get("parent_id", ""),
                "parent_name": data.get("parent_name", ""),
                "layout_x": f"{x:.12g}", "layout_y": f"{y:.12g}",
            })

    with (output_dir / "budapest_lcc_edges.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        fields = ["edge_id", "source", "target", "source_original",
                  "target_original", "connectivity_weight", "metric_length", "n_cells"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for index, edge in enumerate(retained):
            source, target = edge["source_original"], edge["target_original"]
            writer.writerow({
                "edge_id": index, "source": local_id[source], "target": local_id[target],
                "source_original": source, "target_original": target,
                "connectivity_weight": f"{edge['weight']:.12g}",
                "metric_length": f"{length:.12g}", "n_cells": cells,
            })

    with (output_dir / "budapest_lcc_fem.txt").open("w", encoding="ascii") as stream:
        stream.write(f"{len(ordered)} {len(retained)}\n")
        for node in ordered:
            x, y = coordinates[node]
            stream.write(f"{x:.12g} {y:.12g}\n")
        for edge in retained:
            stream.write(f"{local_id[edge['source_original']]} "
                         f"{local_id[edge['target_original']]} {length:.12g} {cells}\n")

    (output_dir / "budapest_lcc_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    source, processed = summary["source"], summary["processed"]
    weights = processed["weight"]
    lines = [
        "Budapest connectome preprocessing summary", "", "Source graph",
        f"  vertices:             {source['vertices']}",
        f"  edges:                {source['edges']}",
        f"  connected components: {source['connected_components']}",
        f"  isolated vertices:    {source['isolated_vertices']}",
        f"  self-loops:           {source['self_loops']}", "", "Processed graph",
        "  selection:            largest connected component",
        f"  vertices:             {processed['vertices']}",
        f"  edges before cleanup: {processed['edges_before_self_loop_removal']}",
        f"  removed self-loops:   {processed['removed_self_loops']}",
        f"  retained edges:       {processed['edges']}",
        f"  degree min/mean/max:  {processed['degree']['minimum']} / "
        f"{processed['degree']['mean']:.6g} / {processed['degree']['maximum']}",
        f"  weight min/Q1/med/Q3/max: {weights['minimum']:.6g} / "
        f"{weights['first_quartile']:.6g} / {weights['median']:.6g} / "
        f"{weights['third_quartile']:.6g} / {weights['maximum']:.6g}", "",
        "FEM conversion", f"  metric length per edge: {length:.12g}",
        f"  cells per edge:          {cells}",
        "  coordinates:            circular layout (not anatomical)",
        "  connectivity weights:    preserved in budapest_lcc_edges.csv", "",
        "Validation", "  GraphML/CSV topology: matching",
        "  GraphML/CSV weights:  matching",
    ]
    (output_dir / "budapest_lcc_summary.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = arguments()
    if args.edge_length <= 0.0 or args.cells_per_edge < 1:
        raise ValueError("Edge length must be positive and cells at least one")
    nodes, edges = read_graphml(args.graphml)
    csv_rows = read_csv(args.csv)
    validate_and_merge(nodes, edges, csv_rows)
    components, adjacency = connected_components(nodes, edges)
    selected = set(components[0])
    retained = [edge for edge in edges
                if edge["source_original"] in selected
                and edge["target_original"] in selected
                and edge["source_original"] != edge["target_original"]]
    retained.sort(key=lambda edge: (id_key(edge["source_original"]),
                                    id_key(edge["target_original"])))
    summary = make_summary(nodes, edges, components, adjacency, selected, retained)
    write_outputs(args.output_dir, nodes, selected, retained, summary,
                  args.edge_length, args.cells_per_edge)
    print("GraphML and CSV representations match.")
    print(f"Prepared largest connected component: {len(selected)} vertices, "
          f"{len(retained)} edges.")
    print(f"Output written to {args.output_dir}")


if __name__ == "__main__":
    main()
