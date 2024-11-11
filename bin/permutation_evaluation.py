#! /usr/bin/env python

""" Module evaluation (Self-consistency and robustness) """

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import graph_tool.all as gt
import csv
import argparse
import pyintergraph
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

    # read the seed genes:
    original_seeds = []
    for line in open(args.seeds, "r"):
        # lines starting with '#' will be ignored
        if line[0] == "#":
            continue
        # the first column in the line will be interpreted as a seed
        line_data = line.strip().split("\t")
        seed_gene = line_data[0]
        original_seeds.append(seed_gene)

    perturbed_seeds = []
    for k in args.permuted_seeds:
        l_seeds = []
        for line in open(k, "r"):
            # lines starting with '#' will be ignored
            if line[0] == "#":
                continue
            # the first column in the line will be interpreted as a seed
            line_data = line.strip().split("\t")
            seed_gene = line_data[0]
            l_seeds.append(seed_gene)
        perturbed_seeds.append(l_seeds)

    # read the PPI network and create a NetworkX graph
    G_ppi = pyintergraph.gt2nx(gt.load_graph(args.network), labelname="name")

    # G_connected_ppi = G_ppi.subgraph(
    #    max(nx.connected_components(G_ppi), key=len)
    # )

    return (
        reference_candidates,
        original_seeds,
        lists_candidates,
        perturbed_seeds,
        G_ppi,
    )


# =============================================================================


def write_output_tsv_file(data, headers, file_name):

    if len(data[0]) == (len(headers) + 1):
        headers = [""] + headers
    assert len(data[0]) == len(headers)

    with open(file_name, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(data)


# =============================================================================


def retrieval_score_from_removed_gene(
    original_seeds, lists_candidates, perturbed_seeds
):
    """
    Self-consistency measure: Compute the averaged retrieval frequency of missing seed
    genes over all cases of missing seed(s).

    Return:
        scores:                 frequency of seed retrieval for all permutations
        avg_scores:             frequency of seed retrieval averaged over all permutations
        std_scores:             corresponding standard deviation
        scores_normalized:      frequency of seed retrieval for all permutations
                                normalized w.r.t. the size of the perturbed lists
        avg_scores_normalized:  frequency of seed retrieval averaged over all permutations
                                normalized w.r.t. the size of the perturbed lists
        std_scores_normalized:  corresponding standard deviation
    """

    scores = []
    scores_normalized = []

    for i in range(len(perturbed_seeds)):

        removed_genes = [j for j in original_seeds if j not in perturbed_seeds[i]]
        n_removed_genes = len(removed_genes)
        score = 0
        for g in removed_genes:
            if g in lists_candidates[i]:
                score += 1 / n_removed_genes
        scores.append(score)
        candidates = [k for k in lists_candidates[i] if k not in perturbed_seeds[i]]
        if len(candidates) == 0:
            scores_normalized.append(0)
        else:
            scores_normalized.append(score / len(candidates))

    avg_scores = round(np.mean(scores), 4)
    std_scores = round(np.std(scores), 4)

    avg_scores_normalized = round(np.mean(scores_normalized), 4)
    std_scores_normalized = round(np.std(scores_normalized), 4)

    scores_normalized_round = [round(k, 4) for k in scores_normalized]

    return (
        scores,
        avg_scores,
        std_scores,
        scores_normalized_round,
        avg_scores_normalized,
        std_scores_normalized,
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


def topological_measures(reference_candidates, G, lists_candidates):
    """
    Robustness measure: Compute four topological measures to compare the reference module with the
    permutated modules:
        - the size of the largest connected component
        - the number of connected genes in the module
        - the number of edges in the module normalized w.r.t the number of possible connections
        (sum over the degree of all nodes)
        - the modularity of the module (module VS the rest of the network)

    Return:
        l_results: list of results for the 4 topological measures for all permutations
        l_mu:      list of the 4 averaged topological measures (average over permutations)
        l_std:     list of the 4 corresponding standard deviations
        l_zscore:  list of the 4 z-scores for the topological measures
    """

    candidate_network = G.subgraph(reference_candidates)
    nodes_ppi = list(G.nodes())
    n_candidates = len(reference_candidates)
    nodes_ppi_wo_candidate = [n for n in nodes_ppi if n not in reference_candidates]

    # size of the largest connected component (LCC)
    lcc_size = len(max(nx.connected_components(candidate_network), key=len))
    l_random_lcc = []

    # number of connected genes (takes also into account the case of multiple connected components)
    interconnected_genes = n_candidates - len(list(nx.isolates(candidate_network)))
    l_interconnected_genes = []

    # normalized number of interedges
    interedges = candidate_network.number_of_edges()

    n_possible_connections = (
        sum([G.degree(s) for s in reference_candidates]) - interedges
    )
    edgibility = interedges / n_possible_connections
    l_random_edgibility = []

    # modularity (computed as the candidates VS the whole network)
    modularity = nx.community.modularity(
        G, [reference_candidates, nodes_ppi_wo_candidate]
    )
    l_random_modularity = []

    for l in lists_candidates:

        # check if the list is empty
        if not l:
            l_random_lcc.append(0)
            l_interconnected_genes.append(0)
            l_random_edgibility.append(0)
            l_random_modularity.append(0)
            logger.warning("Empty list of permuted candidates")
            continue

        G_sub = nx.subgraph(G, l)
        n_candidates_perturbed = G_sub.number_of_nodes()

        # LCC
        lcc_size_rd = len(max(nx.connected_components(G_sub), key=len))
        l_random_lcc.append(lcc_size_rd)

        # number of connected genes
        interconnected_nodes_rd = n_candidates_perturbed - len(list(nx.isolates(G_sub)))
        l_interconnected_genes.append(interconnected_nodes_rd)

        # normalized number of interedges
        interedges_rd = G_sub.number_of_edges()
        n_possible_connections_rd = sum([G.degree(s) for s in l]) - interedges_rd
        edgibility_rd = interedges_rd / n_possible_connections_rd
        l_random_edgibility.append(round(edgibility_rd, 4))

        # modularity
        nodes_ppi_wo_seeds_rd = [n for n in nodes_ppi if n not in l]
        modularity_rd = nx.community.modularity(G, [l, nodes_ppi_wo_seeds_rd])
        l_random_modularity.append(round(modularity_rd, 4))

    mu_lcc = np.mean(l_random_lcc)
    std_lcc = np.std(l_random_lcc)
    z_lcc = (lcc_size - mu_lcc) / std_lcc

    mu_connected_genes = np.mean(l_interconnected_genes)
    std_connected_genes = np.std(l_interconnected_genes)
    z_connected_genes = (
        interconnected_genes - mu_connected_genes
    ) / std_connected_genes

    mu_interedges = np.mean(l_random_edgibility)
    std_interedges = np.std(l_random_edgibility)
    z_interedges = (edgibility - mu_interedges) / std_interedges

    mu_modularity = np.mean(l_random_modularity)
    std_modularity = np.std(l_random_modularity)
    z_modularity = (modularity - mu_modularity) / std_modularity

    # aggregate results
    l_results = [
        lcc_size,
        interconnected_genes,
        round(edgibility, 4),
        round(modularity, 4),
    ]
    l_results_permuted = [
        l_random_lcc,
        l_interconnected_genes,
        l_random_edgibility,
        l_random_modularity,
    ]
    l_mu = [
        round(mu_lcc, 4),
        round(mu_connected_genes, 4),
        round(mu_interedges, 4),
        round(mu_modularity, 4),
    ]
    l_std = [
        round(std_lcc, 4),
        round(std_connected_genes, 4),
        round(std_interedges, 4),
        round(std_modularity, 4),
    ]
    l_zscore = [
        round(z_lcc, 4),
        round(z_connected_genes, 4),
        round(z_interedges, 4),
        round(z_modularity, 4),
    ]

    return l_results, l_results_permuted, l_mu, l_std, l_zscore


# =============================================================================


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse network files to different formats.",
        epilog="Example: python gt2biopax.py network.gt --namespace entrez",
    )
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
        "--seeds",
        help="The original seed file.",
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
        "--permuted_seeds",
        help="The permuted seed files in gt format.",
        type=str,
        required=True,
        nargs="+",
    )

    parser.add_argument(
        "--network",
        help="The reference network in gt format.",
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

    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    print(args)

    reference_candidates, original_seeds, lists_candidates, perturbed_seeds, G = (
        read_input(args)
    )

    # SELF-CONSISTENCY - RETRIEVAL SCORES
    (
        scores,
        avg_scores,
        std_scores,
        scores_normalized,
        avg_scores_normalized,
        std_scores_normalized,
    ) = retrieval_score_from_removed_gene(
        original_seeds, lists_candidates, perturbed_seeds
    )

    # ROBUSTNESS - JACCARD INDEX
    scores_Jaccard, avg_score_Jaccard, std_score_Jaccard = jaccard_index(
        lists_candidates, reference_candidates
    )

    # ROBUSTNESS - TOPOLOGICAL SIMILARITY OF THE MODULES
    l_results, l_results_permuted, l_mu, l_std, l_zscore = topological_measures(
        reference_candidates, G, lists_candidates
    )

    # CREATE TABLES WITH RESULTS
    data_headers = [
        "Retrieval frequency",
        "Normalized retrieval frequency",
        "Jaccard index",
        "LCC size",
        "Number connected genes",
        "Interedges",
        "Modularity",
    ]

    # create a table with the detailed results of all permutations
    data_full = [
        scores,
        scores_normalized,
        scores_Jaccard,
        l_results_permuted[0],
        l_results_permuted[1],
        l_results_permuted[2],
        l_results_permuted[3],
    ]

    file_name_full = f"{args.prefix}.permutation_evaluation_detailed.tsv"
    write_output_tsv_file(np.transpose(data_full), data_headers, file_name_full)

    # create a table with the summarized results of the permutations
    data_summary = [
        ["average", "standard deviation", "z-score"],
        [avg_scores, std_scores, "-"],
        [avg_scores_normalized, std_scores_normalized, "-"],
        [avg_score_Jaccard, std_score_Jaccard, "-"],
        [l_mu[0], l_std[0], l_zscore[0]],
        [l_mu[1], l_std[1], l_zscore[1]],
        [l_mu[2], l_std[2], l_zscore[2]],
        [l_mu[3], l_std[3], l_zscore[3]],
    ]

    file_name_summary = f"{args.prefix}.permutation_evaluation_summary.tsv"
    write_output_tsv_file(np.transpose(data_summary), data_headers, file_name_summary)

    # write multiqc summary
    with open(f"{args.prefix}.permutation_multiqc_summary.tsv", "w") as f:
        f.write(
            "id\tavg_jaccard_index\trediscovery_rate\tnormalized_rediscovery_rate\n"
        )
        f.write(
            f"{args.prefix}\t{avg_score_Jaccard}\t{avg_scores}\t{avg_scores_normalized}\n"
        )


if __name__ == "__main__":
    sys.exit(main())
