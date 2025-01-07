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


def save_multiqc(g, stem):
    with open("input_network_mqc.tsv", "w") as file:
        file.write("Network\tNodes\tEdges\n")
        file.write(f"{stem}\t{g.num_vertices()}\t{g.num_edges()}\n")


def save_diamond(g, stem):
    with open(f"{stem}.diamond.csv", "w") as file:
        writer = csv.writer(file, lineterminator="\n")
        for e in g.iter_edges():
            writer.writerow(
                [g.vp["name"][e[0]], g.vp["name"][e[1]]]
            )  # raw edge values are hashed vertex names


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
            writer.writerow(
                [g.vp["name"][e[0]], g.vp["name"][e[1]]]
            )  # raw edge values are hashed vertex names


def save_rwr(g, stem):
    with open(f"{stem}.rwr.csv", "w") as file:
        writer = csv.writer(file, lineterminator="\n")
        for e in g.iter_edges():
            writer.writerow(
                [g.vp["name"][e[0]], g.vp["name"][e[1]]]
            )  # raw edge values are hashed vertex names


def save_sca(g, stem):
    with open(f"{stem}.sca.tsv", w) as file:
        writer = csv.writer(file, delimiter="/t")
        for e in g.iter_edges():
            writer.writerow(
                [
                    g.vp["name"][e[0]],
                    g.vp["name"][e[0]],
                    g.vp["name"][e[1]],
                    g.vp["name"][e[1]],
                ]
            )


def save(g, stem, format):
    """
    Saves a graph_tools Graph object in a specified format
    """
    if format == "gt":
        save_gt(g=g, stem=stem)
        save_multiqc(g=g, stem=stem)
    elif format == "diamond":
        save_diamond(g=g, stem=stem)
    elif format == "domino":
        save_domino(g=g, stem=stem)
    elif format == "robust":
        save_robust(g=g, stem=stem)
    elif format == "rwr":
        save_rwr(g=g, stem=stem)
    else:
        logger.critical(f"Unknown output format: {format}")
        sys.exit(1)


def load(file_in, extension):
    """
    Loads a graph_tools Graph object.
    """
    if extension in [".gt", ".graphml", ".xml", ".dot", ".gml"]:
        return gt.load_graph(str(file_in))
    else:
        return gt.load_graph_from_csv(str(file_in))


def parse_format(file_in, format):
    stem = Path(file_in).stem
    extension = Path(file_in).suffix
    logger.debug(f"{stem=}")
    logger.debug(f"{extension=}")

    g = load(file_in=file_in, extension=extension)
    logger.debug(f"{g=}")

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
        help="Output format (default gt). If format it gt, a summary file for multiqc will be generated as well.",
        choices=("gt", "diamond", "domino", "robust", "rwr"),
        default="gt",
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
    parse_format(args.file_in, args.format)


if __name__ == "__main__":
    sys.exit(main())
