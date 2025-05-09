#!/usr/bin/env python

"""Randomize an input network file (degree-preserving randomization)."""

import argparse
import logging
import sys
from pathlib import Path
import graph_tool.all as gt
import util


logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse the modules of different tools.",
        epilog="Example: python randomize_network.py --network network.txt",
    )
    parser.add_argument(
        "-g",
        "--network",
        help="Path to the network file used for module generation.",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the output file.",
        type=str,
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
    logger.debug(f"{args=}")

    graph = util.load_graph(str(args.network))

    gt_graph = gt.Graph(graph, prune=True)
    n_failed = gt.random_rewire(
        gt_graph, model="constrained-configuration", n_iter=100, edge_sweep=True
    )

    gt_graph.save(args.output)
    if n_failed > 0:
        logger.warning(
            f"Number of rejected edge moves (due to parallel edges or self-loops): {n_failed}"
        )


if __name__ == "__main__":
    sys.exit(main())
