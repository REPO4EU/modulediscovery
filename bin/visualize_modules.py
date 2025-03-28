#!/usr/bin/env python

"""Takes a gt file and saves it as pdf, png, svg, and html."""


import argparse
import logging
import sys
from pathlib import Path

import util

import graph_tool.all as gt
import pandas as pd
import networkx as nx
import pyintergraph
from pyvis.network import Network

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse the modules of different tools.",
        epilog="Example: python save_modules.py -m module1.gt -p 'module1'",
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
        "-n",
        "--max_nodes",
        help="If the module has more nodes it will not be visualized.",
        type=int,
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

    # check if module is small enough
    if g.num_vertices() > args.max_nodes:
        logger.warn(
            f"Not visualizing {args.prefix} since it has too many ({g.num_vertices()}) nodes"
        )
        sys.exit(0)

    # get vertex properties as dataframe
    vp_df = util.vp2df(g)

    # color the seed genes red
    g.vp["color"] = g.new_vertex_property("string")
    for v in g.vertices():
        if g.vertex_properties["is_seed"][v]:
            g.vp["color"][v] = "red"  # red for seed genes
        else:
            g.vp["color"][v] = "blue"  # blue for added genes

    # calculate the layout
    pos = gt.sfdp_layout(g)

    # save as pdf, png, svg
    for format in ["pdf", "png", "svg"]:
        gt.graph_draw(
            g,
            pos,
            vertex_fill_color=g.vp["color"],
            output_size=(1000, 1000),
            vertex_text=g.vp["name"],
            vorder=g.vp["is_seed"],
            vertex_font_size=12,
            edge_pen_width=3,
            output=f"{args.prefix}.{format}",
        )

    # add x and y properties for consistent html layout
    g.vp["x"] = g.new_vertex_property("double")
    g.vp["y"] = g.new_vertex_property("double")

    for v in g.vertices():
        g.vp["x"][v] = pos[v][0] * 30
        g.vp["y"][v] = pos[v][1] * 30

    # convert to networkx and save as html with pyvis
    nx_graph = pyintergraph.gt2nx(g, labelname="name")

    nt = Network()
    nt.from_nx(nx_graph)

    # add titles
    for node in nt.nodes:
        row = vp_df.loc[node["label"]]
        node[
            "title"
        ] = f"""
            {node['label']}\n
            {'\n'.join(f'{col}: {val}' for col, val in row.items())}
        """

    # turn off physics
    nt.toggle_physics(False)

    # show setting buttons
    nt.show_buttons(filter_=["physics", "interaction", "manipulation"])
    # save as html
    nt.show(f"{args.prefix}.html")


if __name__ == "__main__":
    sys.exit(main())
