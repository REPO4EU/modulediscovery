#!/usr/bin/env python

"""Takes a gt file and saves it as pdf, png, svg, and html."""


import argparse
from collections import defaultdict
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
        "-d",
        "--drugs",
        help="Path to the drug file.",
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
    node_mapping = util.name2index(g)

    # check if module is small enough
    if g.num_vertices() > args.max_nodes:
        logger.warn(
            f"Not visualizing {args.prefix} since it has too many ({g.num_vertices()}) nodes"
        )
        sys.exit(0)

    gene2symbol = {}
    if args.drugs:
        # load the drug file
        drugs = pd.read_csv(args.drugs, sep="\t")
        gene2symbol = {}
        gene2drugs = defaultdict(list)
        drug_set = set(drugs["drug_name"].unique())

        for _, row in drugs.iterrows():
            if row["name"] not in gene2symbol:
                gene2symbol[row["name"]] = row["symbol"]
            if row["name"] not in gene2drugs:
                gene2drugs[row["name"]] = []
            if pd.notna(row["drug_name"]):
                gene2drugs[row["name"]].append(row["drug_name"])
        for name, drug_list in gene2drugs.items():
            if name in node_mapping:
                protein_vertex = node_mapping[name]
            else:
                logger.warn(f"Protein {name} not found in the module")
                continue
            for drug in drug_list:
                if not drug in node_mapping:
                    drug_vertex = g.add_vertex()
                    g.vp["name"][drug_vertex] = drug
                    node_mapping[drug] = drug_vertex
                drug_vertex = node_mapping[drug]
                g.add_edge(protein_vertex, drug_vertex)

    # color the seed genes red
    g.vp["color"] = g.new_vertex_property("string")
    g.vp["label_node"] = g.new_vertex_property("string")
    for v in g.vertices():
        if g.vp["name"][v] in gene2symbol:
            g.vp["label_node"][v] = gene2symbol[g.vp["name"][v]]
        else:
            g.vp["label_node"][v] = g.vp["name"][v]
        if args.drugs and g.vp["name"][v] in drug_set:
            g.vp["color"][v] = "green"  # green for drugs
        else:
            if g.vertex_properties["is_seed"][v]:
                g.vp["color"][v] = "red"  # red for seed genes
            else:
                g.vp["color"][v] = "blue"  # blue for added genes

    # calculate the layout
    pos = gt.sfdp_layout(g)

    # save as pdf, png, svg
    for format in ["pdf", "png", "svg"]:
        size = (3000, 3000) if format == "png" else (1000, 1000)
        if g.num_vertices() > 300:
            font_size = 8 if format == "png" else 2
        elif g.num_vertices() > 200:
            font_size = 12 if format == "png" else 4
        elif g.num_vertices() > 100:
            font_size = 16 if format == "png" else 6
        elif g.num_vertices() > 50:
            font_size = 18 if format == "png" else 8
        elif g.num_vertices() > 25:
            font_size = 24 if format == "png" else 12
        else:
            font_size = 36 if format == "png" else 18
        gt.graph_draw(
            g,
            pos,
            vertex_fill_color=g.vp["color"],
            output_size=size,
            vertex_text=g.vp["label_node"],
            vorder=g.vp["is_seed"],
            vertex_font_size=font_size,
            edge_pen_width=2,
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
    # get vertex properties as dataframe
    vp_df = util.vp2df(g)

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
