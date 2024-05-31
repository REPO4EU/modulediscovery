#!/usr/bin/env python

"""Provide a command line tool to parse the output module of different tools."""


import argparse
import csv
import logging
import sys
import os
import graph_tool.all as gt
from pathlib import Path

logger = logging.getLogger()


def filter_diamond(g, module, filter_column, seeds):
    # Diamond uses a tab separated file format
    g.vp["rank"] = g.new_vertex_property("int")
    g.vp["p_hyper"] = g.new_vertex_property("double")
    with open(module, "r") as file:
        reader = csv.DictReader(file, lineterminator="\n", delimiter="\t")
        for row in reader:
            v = gt.find_vertex(g, g.vp.name, row["DIAMOnD_node"])[0]
            g.vp["rank"][v] = row["#rank"]
            g.vp["p_hyper"][v] = row["p_hyper"]
            g.vp[filter_column][v] = True

        # add seed genes
        for seed in seeds:
            v = gt.find_vertex(g, g.vp.name, seed)
            if v:
                g.vp[filter_column][v[0]] = True
            else:
                logger.warning(
                    f"Did not add seed {seed} since it is not in the network."
                )
    return g


def filter_domino(g, module, filter_column):
    module_ids = []
    with open(module, "r") as file:
        for line in file:
            module_ids += [id.strip("entrez.") for id in line.strip("[]\n").split(", ")]
    gt.map_property_values(
        g.vp.name, g.vp[filter_column], lambda name: name in module_ids
    )
    return g


def filter_robust(g, module, filter_column):
    import numpy as np

    g = gt.load_graph(str(module))
    g.vp.name = g.vp._graphml_vertex_id.copy()
    del g.vp["_graphml_vertex_id"]
    del g.ep["_graphml_edge_id"]
    g.vp[filter_column] = g.new_vertex_property("bool")
    g.vp[filter_column].a = gt.PropertyArray(
        np.ones(len(g), dtype=np.uint8), g.vp[filter_column]
    )
    return g


def filter_rwr(g, module, filter_column):
    g.vp["rank"] = g.new_vertex_property("int")
    g.vp["visiting_probability"] = g.new_vertex_property("double")
    with open(module, "r") as file:
        reader = csv.DictReader(file, lineterminator="\n", delimiter="\t")
        for row in reader:
            v = gt.find_vertex(g, g.vp.name, row["RWR_node"])[0]
            g.vp["rank"][v] = row["#rank"]
            g.vp["visiting_probability"][v] = row["visiting_probability"]
            g.vp[filter_column][v] = True
    return g


def filter_g(g, tool, module, seeds):
    """
    Filters a graph_tools Graph object based on a module of a given tool.
    """
    filter_column = "keep"
    g.vp[filter_column] = g.new_vertex_property("bool")
    if tool == "diamond":
        g = filter_diamond(g, module, filter_column, seeds)
    elif tool == "domino":
        g = filter_domino(g, module, filter_column)
    elif tool == "robust":
        g = filter_robust(g, module, filter_column)
    elif tool == "rwr":
        g = filter_rwr(g, module, filter_column)
    else:
        logger.critical(f"Unknown tool: {tool}")
        sys.exit(1)
    g.set_vertex_filter(g.vp[filter_column])
    g.purge_vertices()
    g.clear_filters()
    del g.vp[filter_column]
    return g


def mark_seeds(g, seeds, property="is_seed"):
    """
    Marks the seed genes in a Graph via a vertex property.
    """
    g.vp[property] = g.new_vertex_property("bool")
    gt.map_property_values(g.vp.name, g.vp[property], lambda name: name in seeds)
    return g


def load(file_in, extension):
    """
    Loads a graph_tools Graph object.
    """
    if extension in [".gt", ".graphml", ".xml", ".dot", ".gml"]:
        return gt.load_graph(str(file_in))
    else:
        return gt.load_graph_from_csv(str(file_in))


def read_seeds(path):
    """
    Loads a set of seeds from a file containing one line per seed gene.
    """
    with open(path, "r") as file:
        seeds = [line.strip() for line in file.readlines() if line.strip()]
    return set(seeds)


def parse_module(file_in, tool, module, seeds_path, output):
    stem = Path(file_in).stem
    extension = Path(file_in).suffix
    logger.debug(f"{stem=}")
    logger.debug(f"{extension=}")

    g = load(file_in=file_in, extension=extension)
    logger.debug(f"{g=}")

    seeds = read_seeds(seeds_path)
    logger.debug(f"{seeds=}")

    g = filter_g(g, tool, module, seeds)
    g = mark_seeds(g, seeds)
    g.save(output)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse the modules of different tools.",
        epilog="Example: python module_parser.py network.gt -t diamond -m module.txt -o module.gt",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Input network.",
    )
    parser.add_argument(
        "-t",
        "--tool",
        help="The tool, that generated the module.",
        choices=("diamond", "domino", "robust", "rwr"),
    )
    parser.add_argument(
        "-m",
        "--module",
        help="Path to the module output.",
        type=Path,
    )
    parser.add_argument(
        "-s",
        "--seeds",
        help="Path to the seeds file used for module generation.",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the parsed output.",
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
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    logger.debug(f"{args=}")
    parse_module(args.file_in, args.tool, args.module, args.seeds, args.output)


if __name__ == "__main__":
    sys.exit(main())
