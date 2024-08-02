#!/usr/bin/env python

"""Provide a command line tool to parse the output module of different tools."""


import argparse
import csv
import logging
import sys
import os
import graph_tool.all as gt
import util
from pathlib import Path
import pandas as pd

logger = logging.getLogger()


def vp2df(g):
    """Convert the vertex properties of a graph to a pandas DataFrame. 'name' will be used as index."""
    # Get all vertex properties
    vertex_props = g.vertex_properties

    # Prepare a dictionary to hold the data for the DataFrame
    data = {
        prop_name: [prop[v] for v in g.vertices()]
        for prop_name, prop in vertex_props.items()
    }

    # Create and return the DataFrame
    df = pd.DataFrame(data)
    df.set_index("name", inplace=True)
    return df


def ep2df(g):
    """Convert the edge properties of a graph to a pandas DataFrame. 'source' and 'target' will be used as index."""
    # Get all edge properties
    edge_props = g.edge_properties

    # Prepare a dictionary to hold the data for the DataFrame
    data = {
        prop_name: [prop[e] for e in g.edges()]
        for prop_name, prop in edge_props.items()
    }
    source_list = []
    target_list = []
    for e in g.iter_edges():
        source_list.append(g.vp["name"][e[0]])
        target_list.append(g.vp["name"][e[1]])

    data["source"] = source_list
    data["target"] = target_list

    # Create and return the DataFrame
    df = pd.DataFrame(data)
    df.set_index(["source", "target"], inplace=True)
    return df


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse the modules of different tools.",
        epilog="Example: python module_parser.py network.gt -t diamond -m module.txt -o module.gt",
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
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.module.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    logger.debug(f"{args=}")

    # load the module file
    g = util.load_graph(str(args.module))

    # save as graphml
    g.save(f"{args.prefix}.graphml")

    # save nodes as tsv
    vp2df(g).to_csv(f"{args.prefix}.nodes.tsv", sep="\t")

    # save edges as tsv
    ep2df(g).to_csv(f"{args.prefix}.edges.tsv", sep="\t")

    # save figures

    vertex_color = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        if g.vertex_properties["is_seed"][v]:
            vertex_color[v] = [1.0, 0.0, 0.0, 1.0]  # red for seed genes
        else:
            vertex_color[v] = [0.0, 0.0, 1.0, 1.0]  # blue for added genes

    # calculate the layout
    pos = gt.sfdp_layout(g)

    # save as pdf
    for format in ["pdf", "png", "svg"]:
        gt.graph_draw(
            g,
            pos,
            vertex_fill_color=vertex_color,
            output_size=(1000, 1000),
            vertex_text=g.vp["name"],
            vorder=g.vp["is_seed"],
            vertex_font_size=12,
            edge_pen_width=3,
            output=f"{args.prefix}.{format}",
        )


if __name__ == "__main__":
    sys.exit(main())
