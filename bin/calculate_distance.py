#!/usr/bin/env python

import graph_tool.all as gt
import pandas as pd
import sys
import numpy as np
import argparse


def calculate_average_distances_all(g):
    all_vertices = list(g.vertices())
    distances = []

    for v in all_vertices:
        if v.out_degree() == 0:
            continue
        dists = gt.shortest_distance(g, source=v, target=all_vertices)
        finite_dists = [d for d in dists if d != 2147483647]
        distances.extend(finite_dists)

    return round(np.mean(distances), 2)


def find_max_shortest_path(graph, gene_file, property_name):
    genes_df = pd.read_csv(gene_file, sep="\t")
    if "gene_id" in genes_df.columns:
        seed_ids = set(genes_df["gene_id"])
    else:
        genes_df = pd.read_csv(gene_file, sep="\t", header=None)
        seed_ids = set(genes_df.iloc[:, 0])
    name_property = graph.vertex_properties[property_name]
    vertex_index = {name_property[v]: v for v in graph.vertices()}

    seeds = [vertex_index[gene_id] for gene_id in seed_ids if gene_id in vertex_index]
    non_seeds = [v for v in graph.vertices() if name_property[v] not in seed_ids]

    max_shortest_path = float("-inf")
    for non_seed in non_seeds:
        shortest_paths = gt.shortest_distance(graph, source=non_seed, target=seeds)
        min_path_length = min(shortest_paths)
        max_shortest_path = max(max_shortest_path, min_path_length)

    return max_shortest_path, seeds


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate Distances")
    parser.add_argument("--module", required=True, help="Input module path")
    parser.add_argument("--nodes", required=True, help="nodes file")
    parser.add_argument("--out", required=True, help="output file")
    parser.add_argument("--id", required=True, help="Id for the output")

    args = parser.parse_args()
    graph_path = args.module
    nodes_file = args.nodes
    out = args.out

    g = gt.load_graph(graph_path)

    average_distance = calculate_average_distances_all(g)
    max_shortest_path, seeds = find_max_shortest_path(g, nodes_file, "name")

    pseudo_diameter, pseudo_diameter_ends = gt.pseudo_diameter(g)

    with open(out, "w") as file:
        file.write(
            "sample\tnodes\tedges\tdiameter\taverage_shortest_path\tmax_dist_to_seed\n"
        )
        file.write(
            f"{args.id}\t{g.num_vertices()}\t{g.num_edges()}\t{pseudo_diameter}\t{average_distance}\t{max_shortest_path}\n"
        )
