#!/usr/bin/env python

from graph_tool.all import *
import sys
import argparse


def gt_to_tsv(input_file, output_file):
    
    graph = load_graph(input_file)

    name_property = graph.vertex_properties["name"]

    with open(output_file, "w") as output_file:
        output_file.write("gene_id\n")

        for vertex in graph.vertices():
            gene_name = name_property[vertex]
            output_file.write(f"{gene_name}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GT file and generate TSV.")
    parser.add_argument("--input", required=True, help="Input GT file path")
    parser.add_argument("--output", required=True, help="Output TSV file path")

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output

    gt_to_tsv(input_file, output_file)
