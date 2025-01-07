#! /usr/bin/env python
"""
# -----------------------------------------------------------------------
# encoding: utf-8

# SeedConnectorAlgorithm.py
# Ruisheng Wang
# Last Modified: 2024-10-25

# This code runs the seed coonector algorithm as described in
#
# Network-Based Disease Module Discovery by a Novel Seed Connector Algorithm with Pathobiological Implications, J Mol Biol. 2018 Sep 14;430(18 Pt A):2939-2950.
# by Rui-Sheng Wang & Joseph Loscalzo
#
# -----------------------------------------------------------------------
"""
import sys

import networkx as nx


def print_usage():

    print(" ")
    print(
        "        usage: ./SCA network_file seed_file outfile_name1 (optional) and outfile_name2 (optional)"
    )
    print("        -----------------------------------------------------------------")
    print(
        "        network_file : The edgelist must be provided as any delimiter-separated"
    )
    print(
        "                       table. The first and third columns are entrez IDs and the"
    )
    print("                        second and the fourth columns are gene symbols")
    print(
        "                       Each row will be interpreted as an interaction gene1 <==> gene2"
    )
    print(
        "        seed_file    : A single column containing the seed genes (if table contains"
    )
    print("                       more than one column they must be tab-separated;")
    print("                       the first column will be used only)")
    print(
        "        outfile_name1 :The predicted gene list will be saved under this file name"
    )
    print(
        '                       by default the outfile_name is set to "Predicted_genelist.txt"'
    )
    print(
        "        outfile_name2 :The edge list of the final module will be saved under this file name"
    )
    print(
        '                       by default the outfile_name is set to "Network_module.txt"'
    )
    print(" ")


def check_input_style(input_list):

    try:
        network_edgelist_file = input_list[1]
        seeds_file = input_list[2]
        # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return

    outfile1 = "Predicted_genelist.txt"
    outfile2 = "Network_module.txt"

    if len(input_list) >= 4:
        outfile1 = input_list[3]
    if len(input_list) == 5:
        outfile2 = input_list[4]

    return network_edgelist_file, seeds_file, outfile1, outfile2


def read_input(network_file, seed_file):

    ## read the network file
    file1 = open(network_file, "r")
    G = nx.Graph()
    for line in file1:
        tmp = line.split()
        G.add_node(tmp[1])
        G.add_node(tmp[3])
        G.add_edge(tmp[1], tmp[3])
    file1.close()

    ## read the seed gene list
    file2 = open(seed_file, "r")
    seed_genes = []
    for line in file2:
        tmp = line.split()
        if G.has_node(tmp[0]):
            seed_genes.append(tmp[0])
    file2.close()

    return G, seed_genes


def output(G, seed_genes, added_nodes, outfile1, outfile2):

    ## output the seeds and the connectors

    fo1 = open(outfile1, "w")
    for i in range(len(seed_genes)):
        fo1.write("%s\t%s\n" % (seed_genes[i], "seed"))

    for i in range(len(added_nodes)):
        fo1.write("%s\t%s\n" % (added_nodes[i], "connector"))
    fo1.close()

    ## output the final module (edge list)
    fo2 = open(outfile2, "wb")
    predicted_genelist = seed_genes + added_nodes
    Module = G.subgraph(predicted_genelist)
    nx.write_edgelist(Module, fo2, delimiter="\t")
    fo2.close()

    return 1


def SCA(G, seed_genes):

    added_nodes = []
    temp_seed_pool = seed_genes
    candidate_genes = []

    for i in range(len(temp_seed_pool)):
        N = G.neighbors(temp_seed_pool[i])
        NB = list(N)
        L = len(NB)
        for j in range(L):
            if NB[j] not in candidate_genes and NB[j] not in temp_seed_pool:
                candidate_genes.append(NB[j])
    while 1:
        SubNet = G.subgraph(temp_seed_pool)
        Gcc = sorted(nx.connected_components(SubNet), key=len, reverse=True)
        H = G.subgraph(Gcc[0])
        LCC = len(set(H.nodes()).intersection(seed_genes))
        candidate_ranks = [0 for i in range(len(candidate_genes))]
        for k in range(len(candidate_genes)):
            TmpGene = []
            TmpGene = temp_seed_pool[:]
            TmpGene.append(candidate_genes[k])
            TmpNet = G.subgraph(TmpGene)
            Gcc = sorted(nx.connected_components(TmpNet), key=len, reverse=True)
            H1 = G.subgraph(Gcc[0])
            LCC1 = len(set(H1.nodes()).intersection(seed_genes))
            candidate_ranks[k] = LCC1 - LCC
            del TmpGene[:]
        idx = [s for s, x in enumerate(candidate_ranks) if x == max(candidate_ranks)]
        if max(candidate_ranks) == 0:
            break

        s = -1
        r = -1
        for t in range(len(idx)):
            L1 = len(list(G.neighbors(candidate_genes[idx[t]])))
            ratio = len(
                set(temp_seed_pool).intersection(G.neighbors(candidate_genes[idx[t]]))
            ) / float(L)
            if r < ratio:
                r = ratio
                s = t
        temp_seed_pool.append(candidate_genes[idx[s]])
        added_nodes.append(candidate_genes[idx[s]])
        print(candidate_genes[idx[s]])
        candidate_genes.remove(candidate_genes[idx[s]])

    return added_nodes


def main():

    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] file providing the network in the form of an edgelist
    #     (tab-separated table, columns 2 & 4 will be used; Columns 1 & 3 are entrez IDs)
    #
    # [2] file with the seed genes (if table contains more than one
    #     column they must be tab-separated; the first column will be
    #     used only)
    #
    # [3] output file for seed genes and predicted seed connectors
    #
    # [4] output file for the final module (edge list)

    # check if input style is correct
    input_list = sys.argv
    network_edgelist_file, seeds_file, outfile1, outfile2 = check_input_style(
        input_list
    )
    G, seed_genes = read_input(network_edgelist_file, seeds_file)
    added_nodes = SCA(G, seed_genes)
    output(G, seed_genes, added_nodes, outfile1, outfile2)

    print("\n The results have been saved to '%s' and '%s' \n" % (outfile1, outfile2))


if __name__ == "__main__":

    main()
