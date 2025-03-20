#!/usr/bin/env python

import graph_tool.all as gt
import sys
import argparse
import logging

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(description="Topology analysis")
    parser.add_argument("--module", required=True, help="Input module path")
    parser.add_argument("--out", required=True, help="Output file")
    parser.add_argument("--id", required=True, help="Id for the output")
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args()
    graph_path = args.module
    out = args.out

    g = gt.load_graph(graph_path)

    # calculate number of seeds and max distance to closest seed
    if "is_seed" in g.vp:
        seeds = []
        added_nodes = []
        for v in g.vertices():
            if g.vp["is_seed"][v] == 1:
                seeds.append(v)
            else:
                added_nodes.append(v)

        num_seeds = len(seeds)

        max_dist_to_seed = 0
        for added_node in added_nodes:
            shortest_paths = gt.shortest_distance(g, source=added_node, target=seeds)
            min_path_length = min(shortest_paths)
            max_dist_to_seed = max(max_dist_to_seed, min_path_length)
    else:
        num_seeds = ""
        max_dist_to_seed = ""

    # connected components
    component_labels, component_sizes = gt.label_components(g)
    num_components = len(component_sizes)
    largest_component = max(component_sizes)

    # diameter
    pseudo_diameter, pseudo_diameter_ends = gt.pseudo_diameter(g)

    # number of isolated nodes
    num_isolated_nodes = len([v for v in g.vertices() if g.vertex(v).out_degree() == 0])

    # write output
    with open(out, "w") as file:
        file.write(
            "\t".join(
                [
                    "sample",
                    "nodes",
                    "edges",
                    "seeds",
                    "max_dist_to_seed",
                    "diameter",
                    "components",
                    "largest_component",
                    "isolated_nodes",
                ]
            )
            + "\n"
        )
        file.write(
            f"{args.id}\t"
            f"{g.num_vertices()}\t"
            f"{g.num_edges()}\t"
            f"{num_seeds}\t"
            f"{max_dist_to_seed}\t"
            f"{pseudo_diameter}\t"
            f"{num_components}\t"
            f"{largest_component}\t"
            f"{num_isolated_nodes}\n"
        )


if __name__ == "__main__":
    sys.exit(main())
