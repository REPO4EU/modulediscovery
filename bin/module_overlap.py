#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate overlap matrices",
        epilog="Example: python module_overlap.py --ids id1 id2 --inputs file1 file2",
    )
    parser.add_argument(
        "--ids",
        type=str,
        help="IDs to name the columns/rows of the output matrix",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--inputs",
        help="Input files containing node lists (one line per row). Number has to match the number of ids.",
        type=Path,
        nargs="+",
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


def read_node_set(file, skip_header=True):
    """Read nodes from a file."""
    with open(file, "r") as f:
        if skip_header:
            next(f)
        return set([line.strip() for line in f.readlines()])


def pairwise_matrix(sets, set_names, function):
    df = pd.DataFrame(columns=set_names)
    for set_1, name_1 in zip(sets, set_names):
        row = {name_2: function(set_1, set_2) for set_2, name_2 in zip(sets, set_names)}
        df.loc[len(df)] = row
    df.insert(1, "ID", set_names, True)
    df.set_index("ID", inplace=True)
    return df


def jaccard(set_1, set_2):
    return len(set_1.intersection(set_2)) / len(set_1.union(set_2))


def shared_nodes(set_1, set_2):
    return len(set_1.intersection(set_2))


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    assert len(args.ids) == len(args.inputs)

    module_map = {}
    for id, file in zip(args.ids, args.inputs):
        module_map[id] = read_node_set(file)
    module_map = dict(sorted(module_map.items()))
    ids, modules = zip(*module_map.items())

    jaccard_df = pairwise_matrix(modules, ids, jaccard)
    jaccard_df.to_csv("jaccard_similarity_matrix_mqc.tsv", sep="\t")

    shared_nodes_df = pairwise_matrix(modules, ids, shared_nodes)
    shared_nodes_df.to_csv("shared_nodes_matrix_mqc.tsv", sep="\t")


if __name__ == "__main__":
    sys.exit(main())
