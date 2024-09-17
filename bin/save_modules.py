#!/usr/bin/env python

"""Takes a gt file and saves it as graphml and tsv."""


import argparse
import logging
import sys
from pathlib import Path

import util

import graph_tool.all as gt
import pandas as pd

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
    vp_df = util.vp2df(g)
    vp_df.to_csv(f"{args.prefix}.nodes.tsv", sep="\t")

    # save edges as tsv
    util.ep2df(g).to_csv(f"{args.prefix}.edges.tsv", sep="\t")


if __name__ == "__main__":
    sys.exit(main())
