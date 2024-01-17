#!/usr/bin/env python

"""Provide a command line tool to parse different network file formats."""


import argparse
import csv
import logging
import sys
import os
import graph_tool.all as gt
from pathlib import Path

logger = logging.getLogger()


def save_gt(g, stem):
    g.save(f"{stem}.gt")


def save_diamond(g, stem):
    with open(f"{stem}.diamond.csv", "w") as file:
        writer = csv.writer(file, lineterminator="\n")
        for e in g.iter_edges():
            writer.writerow([g.vp["name"][e[0]], g.vp["name"][e[1]]])  # raw edge values are hashed vertex names


def save_domino(g, stem):
    with open(f"{stem}.domino.sif", "w") as file:
        writer = csv.writer(file, lineterminator="\n", delimiter="\t")
        writer.writerow(["node_1", "type", "node_2"])  # write header
        for e in g.iter_edges():
            writer.writerow(
                [
                    f"entrez.{g.vp['name'][e[0]]}",
                    "ppi",
                    f"entrez.{g.vp['name'][e[1]]}",
                ]
            )  # raw edge values are hashed vertex names


def save_robust(g, stem):
    with open(f"{stem}.robust.tsv", "w") as file:
        writer = csv.writer(file, lineterminator="\n", delimiter="\t")
        for e in g.iter_edges():
            writer.writerow([g.vp["name"][e[0]], g.vp["name"][e[1]]])  # raw edge values are hashed vertex names


def save(g, stem, format):
    """
    Saves a graph_tools Graph object in a specified format
    """
    if format == "gt":
        save_gt(g=g, stem=stem)
    elif format == "diamond":
        save_diamond(g=g, stem=stem)
    elif format == "domino":
        save_domino(g=g, stem=stem)
    elif format == "robust":
        save_robust(g=g, stem=stem)
    else:
        logger.critical(f"Unknown output format: {format}")
        sys.exit(1)


def filter_diamond(g, module, filter_column):
    # Diamond uses a tab separated file format
    g.vp["diamond_rank"] = g.new_vertex_property("int")
    g.vp["diamond_p_hyper"] = g.new_vertex_property("double")
    with open(module, "r") as file:
        reader = csv.DictReader(file, lineterminator="\n", delimiter="\t")
        for row in reader:
            v = gt.find_vertex(g, g.vp.name, row["DIAMOnD_node"])[0]
            g.vp["diamond_rank"][v] = row["#rank"]
            g.vp["diamond_p_hyper"][v] = row["p_hyper"]
            g.vp[filter_column][v] = True
    return g


def filter_domino(g, module, filter_column):
    module_ids = []
    with open(module, "r") as file:
        for line in file:
            module_ids += [id.strip("entrez.") for id in line.strip("[]\n").split(", ")]
    gt.map_property_values(g.vp.name, g.vp[filter_column], lambda name: name in module_ids)
    return g


def filter_robust(g, module, filter_column):
    import numpy as np

    g = gt.load_graph(str(module))
    g.vp.name = g.vp._graphml_vertex_id.copy()
    del g.vp["_graphml_vertex_id"]
    del g.ep["_graphml_edge_id"]
    g.vp[filter_column] = g.new_vertex_property("bool")
    g.vp[filter_column].a = gt.PropertyArray(np.ones(len(g), dtype=np.uint8), g.vp[filter_column])
    return g


def filter_g(g, format, module):
    """
    Filters a graph_tools Graph object based on a module file in a given format
    """
    filter_column = "keep"
    g.vp[filter_column] = g.new_vertex_property("bool")
    if format == "gt":
        logger.warning("Module file given, but format is 'gt'. Network was not filtered based on module...")
        return g
    elif format == "diamond":
        g = filter_diamond(g, module, filter_column)
    elif format == "domino":
        g = filter_domino(g, module, filter_column)
    elif format == "robust":
        g = filter_robust(g, module, filter_column)
    else:
        logger.critical(f"Unknown output format: {format}")
        sys.exit(1)
    g.set_vertex_filter(g.vp[filter_column])
    g.purge_vertices()
    g.clear_filters()
    del g.vp[filter_column]
    return g


def load(file_in, extension):
    """
    Loads a graph_tools Graph object.
    """
    if extension in [".gt", ".graphml", ".xml", ".dot", ".gml"]:
        return gt.load_graph(str(file_in))
    else:
        return gt.load_graph_from_csv(str(file_in))


def parse_format(file_in, format, module=None):
    stem = Path(file_in).stem
    extension = Path(file_in).suffix
    logger.debug(f"{stem=}")
    logger.debug(f"{extension=}")

    g = load(file_in=file_in, extension=extension)
    logger.debug(f"{g=}")
    if module:
        g = filter_g(g, format, module)
        stem = format
        format = "gt"
    save(g=g, stem=stem, format=format)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse network files to different formats.",
        epilog="Example: python graph_tools.py network.csv -f gt",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Input network.",
    )
    parser.add_argument(
        "-f",
        "--format",
        help="Output format (default gt).",
        choices=("gt", "diamond", "domino", "robust"),
        default="gt",
    )
    parser.add_argument(
        "-m",
        "--module",
        help="Path to the module output. If this is given, the output will be gt. The -f flag is still used to determine the module format. Network must be in .gt format.",
        type=Path,
        required=False,
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
    parse_format(args.file_in, args.format, args.module)


if __name__ == "__main__":
    sys.exit(main())
