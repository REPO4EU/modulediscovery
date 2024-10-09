#!/usr/bin/env python

import os
import graph_tool.all as gt
from pathlib import Path
import sys
import json
import requests
import argparse
import csv


def load_nodes(graph):
    ##returns a list with the name of the nodes in the JSON format needed for drugstone
    nodes = []
    name = graph.vp["name"]
    seeds = graph.vp["is_seed"]
    for id, is_seed in zip(name, seeds):
        if is_seed:
            group = "seedNode"
        else:
            group = "default"
        nodes.append({"id": id, "group": group})
    return json.dumps(nodes)


def load_edges(graph):
    ##returns a list of edges in the JSON format needed for drugstone
    edges = []
    for source, target in graph.iter_edges():
        nodes = graph.vp["name"]
        source_node = nodes[source]
        target_node = nodes[target]
        edges.append({"from": source_node, "to": target_node})
    return json.dumps(edges)


def send_requests(nodes, edges, identifier):
    ## creates and sends the request to the drugstone api

    url = "https://api.drugst.one/create_network"
    data = {
        "network": {"nodes": nodes, "edges": edges},
        "config": {"identifier": identifier, "autofillEdges": False},
    }
    post_request = requests.post(url, json=data)
    id = post_request.text.strip('"')
    result_url = f"https://drugst.one?id={id}"
    get_r = requests.get(result_url)
    return get_r.url


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="generates a drugstone link for the module",
        epilog="Example: python drugstone.py -m module.gt -o module.drugstonelink.txt -i Entrez",
    )
    parser.add_argument(
        "-m",
        "--module",
        help="Path to the module output.",
        type=Path,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix to name the output files.",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--id_space",
        required=True,
        help="Input Identifier",
        type=str,
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    file = args.module
    id_space = args.id_space
    graph = gt.load_graph(str(file))
    nodes = load_nodes(graph)
    edges = load_edges(graph)
    link = send_requests(nodes, edges, id_space)
    header = ["Module_id", "drugstone_link", "link_raw"]
    data = [f"{args.prefix}", f"<a href={link}>{args.prefix}</a>", f"{link}"]
    with open(f"{args.prefix}.drugstonelink.tsv", "w") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(header)
        writer.writerow(data)


if __name__ == "__main__":
    main()
