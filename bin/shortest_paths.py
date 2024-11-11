#! /usr/bin/env python

import logging
import os
import pickle
import sys
import networkx as nx


def get_shortest_paths(graph, dump_file):
    if not os.path.exists(dump_file):
        lengths = dict(nx.shortest_path_length(graph))
        if dump_file is not None:
            logging.warning("Pickled file not found: {}".format(dump_file))
            pickle.dump(lengths, open(dump_file, "wb"))
    else:
        logging.info("Using existing dump file: {}".format(dump_file))


def parse_network(network_file):  # , id_mapping_file=None
    network = nx.read_edgelist(network_file, delimiter=",")
    component_nodes = max(nx.connected_components(network), key=len)
    network = nx.subgraph(network, component_nodes)
    # if id_mapping_file is not None:
    #     with open(id_mapping_file) as mapping:
    #         mapping = {k: gene for k, gene in (l.strip().split("\t") for l in mapping)}
    #     network = nx.relabel_nodes(network, mapping)
    return network


def main():
    logging.basicConfig(filename="shortest_paths.log", level=logging.DEBUG)  # INFO

    if len(sys.argv) < 2:
        print("Usage: python shortest_paths.py <network_data> [<dump_file>]")
        sys.exit(1)

    network_file = sys.argv[1]
    network = parse_network(network_file)
    dump_file = sys.argv[2] if len(sys.argv) > 2 else "shortest_paths.pkl"
    get_shortest_paths(network, dump_file)


if __name__ == "__main__":
    main()
