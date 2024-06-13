#!/usr/bin/env python

"""Provide a command line tool for pipeline input validation."""


import argparse
import csv
import logging
import sys
import os
import graph_tool.all as gt
from pathlib import Path
import util

logger = logging.getLogger()


def save_gt(g, stem):
    g.save(f"{stem}.gt")


def load_graph(file):
    stem = Path(file).stem
    extension = Path(file).suffix
    logger.debug(f"{stem=}")
    logger.debug(f"{extension=}")
    if extension in [".gt", ".graphml", ".xml", ".dot", ".gml"]:
        return gt.load_graph(str(file))
    else:
        return gt.load_graph_from_csv(str(file))


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Check seeds.",
        epilog="Example: python input_check.py --network network.csv --seeds seeds.txt",
    )
    parser.add_argument(
        "-n",
        "--network",
        type=Path,
        help="The input network.",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix for naming the output files.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s",
        "--seeds",
        help="Path to the seeds file used for module generation.",
        type=Path,
        required=True,
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
    if not args.network.is_file():
        logger.error(f"The given input file {args.network} was not found!")
        sys.exit(2)
    if not args.seeds.is_file():
        logger.error(f"The given input file {args.network} was not found!")
        sys.exit(2)
    logger.debug(f"{args=}")

    g = util.load_graph(str(args.network))
    seeds = util.read_seeds(str(args.seeds))
    logger.debug(f"{seeds=}")

    seeds_keep = []
    seeds_remove = []

    name2index = util.name2index(g)
    for seed in seeds:
        if seed in name2index:
            seeds_keep.append(seed)
        else:
            seeds_remove.append(seed)

    if seeds_keep:
        with open(f"{args.prefix}.tsv", "w") as file:
            for seed in seeds_keep:
                file.write(f"{seed}\n")

    if seeds_remove:
        logger.warning(f"Seeds not in graph: {seeds_remove}")

        with open(f"{args.prefix}.removed.tsv", "w") as file:
            for seed in seeds_remove:
                file.write(f"{seed}\n")

        with open(f"{args.prefix}.multiqc.tsv", "w") as file:
            file.write("Seed file\tSeeds removed\tSeeds remaining\n")
            file.write(f"{args.prefix}\t{len(seeds_remove)}\t{len(seeds_keep)}\n")


if __name__ == "__main__":
    sys.exit(main())
