#! /usr/bin/env python

""" Module evaluation (Self-consistency and robustness) """

import numpy as np
import sys
import graph_tool.all as gt
import argparse
import logging

logger = logging.getLogger()

# =============================================================================


def read_input(args):

    # read the modules from gt format to arrays of genes
    g_init = gt.load_graph(args.module)
    reference_candidates = [g_init.vp["name"][v] for v in g_init.iter_vertices()]

    lists_candidates = []
    for g in args.permuted_modules:
        graph = gt.load_graph(g)
        l_cand = [graph.vp["name"][v] for v in graph.iter_vertices()]
        lists_candidates.append(l_cand)

    return (
        reference_candidates,
        lists_candidates,
    )


# =============================================================================


def jaccard_index(lists_candidates, reference_candidates):
    """
    Robustness measure: Compute the Jaccard index between the reference module and all
    permuted modules.

    Return:
        scores_Jaccard:     Jaccard index between the reference module and all permuted modules separately
        avg_score_Jaccard:  average Jaccard index between the reference module and all permuted modules
        std_score_Jaccard:  standard deviation of Jaccard index between the reference module and all permuted modules
    """

    scores_Jaccard = []
    for list in lists_candidates:
        overlap = len(set(list) & set(reference_candidates))
        scores_Jaccard.append(
            overlap / (len(list) + len(reference_candidates) - overlap)
        )

    avg_score_Jaccard = round(np.mean(scores_Jaccard), 4)
    std_score_Jaccard = round(np.std(scores_Jaccard), 4)
    scores_Jaccard_round = [round(k, 4) for k in scores_Jaccard]

    return scores_Jaccard_round, avg_score_Jaccard, std_score_Jaccard


# =============================================================================


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix to name the output files.",
        type=str,
    )

    parser.add_argument(
        "--module",
        help="The original module file in gt format.",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--permuted_modules",
        help="The permuted module files in gt format.",
        type=str,
        required=True,
        nargs="+",
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

    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    print(args)

    reference_candidates, lists_candidates = read_input(args)

    # ROBUSTNESS - JACCARD INDEX
    scores_Jaccard, avg_score_Jaccard, std_score_Jaccard = jaccard_index(
        lists_candidates, reference_candidates
    )

    # write multiqc summary
    with open(f"{args.prefix}.network_permutation_multiqc_summary.tsv", "w") as f:
        f.write("id\tavg_jaccard_index\n")
        f.write(f"{args.prefix}\t{avg_score_Jaccard}\n")

    # write multiqc jaccard indices
    with open(f"{args.prefix}.network_permutation_multiqc_jaccard.txt", "w") as f:
        f.write(f"{args.prefix}: {scores_Jaccard}\n")


if __name__ == "__main__":
    sys.exit(main())
