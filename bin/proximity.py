#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Emre Guney - emre@stalicla.com - 30/09/2020
#
# ExPERT: Expanded Pathway Perturbation based Repositioning
# of Target Endophenotypes
#
# Implementation of closeness to gene sets from a set of nodes
# in the protein interaction network.
###############################################################
import logging
import os
import pickle
import random
import time
import sys
import numpy
import configparser
import networkx as nx
import pandas as pd
import graph_tool.all as gt
import pyintergraph


def run_proximity(
    drug_to_targets,
    phenotype_to_genes,
    network,
    output_file="./proximity.txt",
    n_random=1000,
    min_bin_size=100,
    seed=51234,
    degree_aware=True,
    phenotype_to_info=None,
    network_dump_file=None,
    skip_no_targets_in_module=False,
):
    """
    Run proximity on each gene set using the provided targets, output is saved
    in a text file
    """
    # Get network
    nodes_network = set(network.nodes())
    # Get shortest paths
    # Get degree binning

    lengths = get_shortest_paths(network_dump_file)

    bins = get_degree_binning(network, min_bin_size, lengths)

    k = 1
    f = init_outFile(output_file, phenotype_to_info)
    for drug, targets in drug_to_targets.items():
        targets = set(targets) & nodes_network

        if len(targets) == 0:
            continue

        for phenotype, genes in phenotype_to_genes.items():
            logging.info("{}/{} {}".format(k, len(phenotype_to_genes), phenotype))
            t1 = time.perf_counter()
            k += 1

            genes = list(sorted(set(genes) & nodes_network))
            if len(genes) == 0:
                # logging.info("Skipping: {}".format(phenotype))
                continue

            d = calculate_min_avg_distance(network, targets, genes, lengths)
            # print(phenotype, len(genes), d)
            random_gene_sets = pick_random_nodes_matching_selected(
                network,
                bins,
                genes,
                n_random=n_random,
                degree_aware=degree_aware,
                seed=seed,
            )
            # print(random_gene_sets[0]) # for testing
            values = numpy.empty(n_random)
            for i, nodes_random in enumerate(random_gene_sets):
                values[i] = calculate_min_avg_distance(
                    network, targets, nodes_random, lengths
                )
            m, s = numpy.mean(values), numpy.std(values)
            if s == 0:
                z = 0.0
            else:
                z = (d - m) / s

            if phenotype_to_info is not None:
                f.write(
                    f"{drug}\t{phenotype}\t{z}\t{phenotype_to_info[0]}\t{phenotype_to_info[1]}\t{d}\t{m}\t{s}\n"
                )
            else:
                f.write(f"{drug}\t{phenotype}\t{z}\t{d}\t{m}\t{s}\n")
            # print(phenotype, len(genes), z, d, m, s)

            f.flush()
            t2 = time.perf_counter()
            logging.info(f"Completed in {t2-t1:.4f}s")

    f.close()
    return None


def get_shortest_paths(dump_file):
    if os.path.exists(dump_file):
        lengths = pickle.load(open(dump_file, "rb"))
    return lengths


def calculate_min_avg_distance(network, nodes_from, nodes_to, lengths):
    """
    Helper function to calculate avg distance to the closest node
    """
    values_outer = []
    for node_from in nodes_from:
        values = []
        vals = lengths[node_from]
        for node_to in nodes_to:
            val = vals[node_to]
            values.append(val)
        d = min(values)
        values_outer.append(d)
    d = numpy.mean(values_outer)
    return d


def get_degree_binning(g, bin_size, lengths=None):
    """
    Helper function to bin nodes based on their degree
    """
    degree_to_nodes = {}
    for node, degree in g.degree():
        # if lengths is given, it will only use those nodes
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)

    values = list(degree_to_nodes)
    values.sort()
    # observed node degrees

    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        # print(i, low, high, len(val))
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins


def pick_random_nodes_matching_selected(
    network,
    bins,
    nodes_selected,
    n_random,
    degree_aware=True,
    connected=False,
    seed=None,
):
    """
    Function to pick random nodes matching the degrees of given nodes.
    bins variable is generated using get_degree_binning (to get bins)
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = list(network.nodes())
    for _ in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(
                nodes_selected, bins, network
            )
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                # nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20):  # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while k < len(nodes_selected):
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(
                        [x for x in network.neighbors(node_random)]
                    )
                    if not node_selected in nodes_random:
                        nodes_random.append(node_selected)
                        k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values


def get_degree_equivalents(seeds, bins, g):
    """
    Order of entries are important due to the removal process
    For reproducibility make sure to feed it with the same ordered list
    """
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes


def init_outFile(output_file, phenotype_to_info):
    f = open(output_file, "w")
    if phenotype_to_info is not None:
        f.write("drug\tphenotype\tz\tmoa\tconsistency\td\tm\ts\n")
    else:
        f.write("drug\tphenotype\tz\td\tm\ts\n")
    return f


def parse_input_file(f, source_column, target_column, prefix):
    source_to_targets = {}
    if source_column:
        data = pd.read_csv(f, sep="\t").set_index(source_column)  # sep was "\t"

        for i, g in data.groupby(source_column):
            tmp = {i: [str(x) for x in g[target_column]]}
            source_to_targets |= tmp
    else:
        data = pd.read_csv(f, sep="\t")
        source_to_targets = {prefix: [str(x) for x in data[target_column]]}

    return source_to_targets


def parse_network(network_file, id_mapping_file=None):
    g = gt.load_graph(network_file)
    network = pyintergraph.gt2nx(g, labelname="name")
    component_nodes = max(nx.connected_components(network), key=len)
    network = nx.subgraph(network, component_nodes)

    if id_mapping_file is not None:
        with open(id_mapping_file) as mapping:
            mapping = {k: gene for k, gene in (l.strip().split("\t") for l in mapping)}
        network = nx.relabel_nodes(network, mapping)
    return network


def parseConfigFile(config_file):
    """REQUIRED FILEDS
    drug_to_target
    drug_column
    target_column

    phenotype_to_gene
    phenotype_column
    gene_column

    network_file
    shortest_path_file = config["network_file"]
    output_file = config["output_file"]

    id_mapping_file = config["PROXIMITY"]["id_mapping_file"]

    OPTIONAL
    n_random: default 1000
    min_bin_size: default 100
    random_seed: default 51234
    shortest_path_file: None
    """

    config = configparser.ConfigParser()
    config.read_dict(
        {
            "PROXIMITY": {
                "shortest_path_file": "None",
                "id_mapping_file": "None",
                "n_random": "1000",
                "min_bin_size": "100",
                "random_seed": "51234",
                "degree_aware": "True",
            }
        }
    )

    config.read(config_file)
    config = config["PROXIMITY"]

    config = {x: config[x] for x in config}
    if config["id_mapping_file"] == "None":
        config["id_mapping_file"] = None
    if config["phenotype_column"] == "None":
        config["phenotype_column"] = None
    if config["degree_aware"] == "True":
        config["degree_aware"] = True
    else:
        config["degree_aware"] = False
    return config


def main():
    logging.basicConfig(filename="proximity.log", level=logging.DEBUG)  # INFO
    config_file = sys.argv[1]
    config = parseConfigFile(config_file)
    drug_to_target = parse_input_file(
        config["drug_to_target"], config["drug_column"], config["target_column"], None
    )

    phenotype_to_genes = parse_input_file(
        config["phenotype_to_gene"],
        config["phenotype_column"],
        config["gene_column"],
        config["prefix"],
    )

    network = parse_network(config["network_file"], config["id_mapping_file"])

    logging.info(
        f"Config file: {config_file}\ngene_set_file: {config['phenotype_to_gene']}\nTargets: {config['drug_to_target']}\n"
    )
    logging.debug(f"{len(phenotype_to_genes)} , {list(phenotype_to_genes.items())[:4]}")
    logging.debug(
        f"{len(network.nodes())} {len(network.edges())} {list(network.edges())[:4]}"
    )

    run_proximity(
        drug_to_target,
        phenotype_to_genes,
        network,
        output_file=config["output_file"],
        n_random=int(config["n_random"]),
        min_bin_size=int(config["min_bin_size"]),
        seed=int(config["random_seed"]),
        degree_aware=config["degree_aware"],
        network_dump_file=config["shortest_paths"],
    )


if __name__ == "__main__":
    main()
