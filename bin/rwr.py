#! /usr/bin/env python


import networkx as nx
import numpy as np
import scipy.sparse as sp
import sys
import csv


# =============================================================================


def print_usage():
    print("")
    print("        usage: ./rwr.py network_edgelist_file seeds_file scaling symmetrical r(optional)")
    print("        -----------------------------------------------------------------")
    print("        network_file : The edgelist must be provided as any delimiter-separated")
    print("                       table. Make sure the delimiter does not exit in gene IDs")
    print("                       and is consistent across the file.")
    print("                       The first two columns of the table will be")
    print("                       interpreted as an interaction gene1 <==> gene2")
    print("        seed_file    : Table containing the seed genes (if table contains")
    print("                       more than one column they must be tab-separated;")
    print("                       the first column will be used only).")
    print("        scaling      : (0 or 1) value to allow (1) the scaling the nodes' visiting")
    print("                       probabilities by the sqrt of their degree or not (0).")
    print("        symmetrical  : (0 or 1) value for the use of the symmetric Markov matrix (1)")
    print("                       (instead of the column-wise normalized Markov matrix (0))")
    print("        r            : damping factor/restart probability for the random walk,")
    print("                       default is set to 0.8.")
    print(" ")


# =============================================================================


def check_input_style(input_list):
    try:
        network_edgelist_file = input_list[1]
        seeds_file = input_list[2]
        scaling = int(input_list[3])
        symmetrical = int(input_list[4])
    # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return

    r = 0.8

    if len(input_list) == 6:
        try:
            r = float(input_list[5])
        except:
            print_usage()
            sys.exit(0)
            return

    sc = ["no_scaling", "scaling"]
    sy = ["columwise", "symmetrical"]

    outfile_name = f"connected_module_rwr_{sc[scaling]}_{sy[symmetrical]}_{r}.txt"

    return network_edgelist_file, seeds_file, scaling, symmetrical, r, outfile_name


# =============================================================================


def read_input(network_file, seed_file):
    """
    Reads the network and the list of seed genes from external files.

    * The edgelist must be provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2

    * The seed genes mus be provided as a table. If the table has more
    than one column, they must be tab-separated. The first column will
    be used only.

    * Lines that start with '#' will be ignored in both cases
    """

    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file, "r"):
        if line[0] == "#":
            continue
        else:
            dialect = sniffer.sniff(line)
            line_delimiter = dialect.delimiter
            break
    if line_delimiter == None:
        print("network_file format not correct")
        sys.exit(0)

    # read the network:
    G = nx.Graph()
    for line in open(network_file, "r"):
        # lines starting with '#' will be ignored
        if line[0] == "#":
            continue
        # The first two columns in the line will be interpreted as an
        # interaction gene1 <=> gene2
        # line_data   = line.strip().split('\t')
        line_data = line.strip().split(line_delimiter)
        node1 = line_data[0]
        node2 = line_data[1]
        G.add_edge(node1, node2)

    G_connected = G.subgraph(max(nx.connected_components(G), key=len))  # extract lcc graph

    # read the seed genes:
    seed_genes = set()
    for line in open(seed_file, "r"):
        # lines starting with '#' will be ignored
        if line[0] == "#":
            continue
        # the first column in the line will be interpreted as a seed
        # gene:
        line_data = line.strip().split("\t")
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)

    return G_connected, seed_genes


# =============================================================================


def colwise_rnd_walk_matrix(G, r, a):
    """
    Compute the Random Walk Matrix (RWM) for a given graph G with teleportation
    probability a and damping factor r using the formula (I-r*M)^-1 where M is
    the column-wise normalized Markov matrix according to M = A D^{-1}

    Parameters:
        G: (networkx graph) input graph
        r: (float) damping factor/restart probability
        a: (float) teleportation probability

    Returns:
        W: (numpy array) RWM of the input graph G

    """

    # get the number of nodes in the graph G
    n = G.number_of_nodes()
    # get the adjacency matrix of graph G
    A = nx.adjacency_matrix(G, sorted(G.nodes()))

    # calculate the first scaling term
    factor = float((1 - a) / n)

    # prepare the second scaling term
    E = np.multiply(factor, np.ones([n, n]))
    A_tele = np.multiply(a, A) + E

    # compute the column-wise normalized Markov matrix
    norm = np.linalg.norm(A_tele, ord=1, axis=0)
    M = np.array(A_tele / norm)

    # mixture of Markov chains
    del A_tele
    del E

    U = np.identity(n, dtype=int)
    H = (1 - r) * M
    H1 = np.subtract(U, H)
    del U
    del M
    del H

    # compute the RWM using the formula (I-r*P)^-1
    W = r * np.linalg.inv(H1)

    return W


# =============================================================================


def symmetric_rnd_walk_matrix(G, r):
    """
    Compute the Random Walk Matrix (RWM) for a given graph G damping factor r
    using the formula (I-r*M_s)^-1 where M_s is the symmetric Markov matrix
    according to M_s = D^{-1/2} A D^{-1/2} = I - Laplace_normalized

    Parameters:
        G: (networkx graph) input graph
        r: (float) damping factor/restart probability

    Returns:
        W: (numpy array) RWM of the input graph G

    """

    # get the number of nodes in the graph G
    n = G.number_of_nodes()

    # generate symmetric Markov matrix
    M_laplace = nx.normalized_laplacian_matrix(G, sorted(G.nodes()))
    Id = sp.identity(n)
    M_s = Id - M_laplace

    H = (1 - r) * M_s
    H1 = np.subtract(Id, H)
    del M_s
    del H

    # compute the RWM using the formula (I-r*M_s)^-1
    W = r * np.linalg.inv(H1.todense())
    W = np.array(W)

    return W


# =============================================================================


def create_mapping_index_entrezID(G):
    """
    Create the dictionaries to map genes' entrez IDs to indices and vice-versa.

    Parameter:
        G: (networkx graph) input graph

    Returns:
        d_entz_idx: dictionary entrez ID to index
        d_idx_entz: dictionare index to entrez ID
    """

    d_idx_entz = {}
    cc = 0
    for entz in sorted(G.nodes()):
        d_idx_entz[cc] = entz
        cc += 1
    d_entz_idx = dict((y, x) for x, y in d_idx_entz.items())

    return d_entz_idx, d_idx_entz


# =============================================================================


def create_scaling_matrix(G):
    """
    Compute the diagonal matrix of the inverse degree of the nodes in graph G.

    Parameter:
        G: (networkx graph) input graph

    Returns:
        Dinvsqrt: diagonal matrix of the inverse degree of the nodes in graph G

    """

    n_nodes = G.number_of_nodes()
    Dinvsqrt = np.zeros([n_nodes, n_nodes])  # basically dimensions of your graph
    cc = 0
    for node in sorted(G.nodes()):
        kn = G.degree(node)
        Dinvsqrt[cc, cc] = np.sqrt(1.0 / kn)
        cc += 1

    return Dinvsqrt


# =============================================================================


def rwr(G, seed_genes, scaling, symmetrical, restart_parameter=0.8, alpha=1.0):
    """
    Perform the random walk process (column-wise or symmetrical, with scaling or not),
    find the visiting probability to each node and determine the top-k ranked genes
    which connect the seed genes.
    The function writes the results (all ranked genes with their visiting probability)
    and outputs the connected disease module.

    Parameters:
        G:                  (networkx graph) input graph
        seed_genes:         (list) seed genes
        scaling:            (boolean) scale the visiting probabilities with the sqrt
                            of the degree of the corresponding node
        symmetrical:        (boolean) compute the symmetric Markov matrix if True
                            the column-wise normalized otherwise
        restart_parameter:  (float) damping factor/restart probability (default value 0.8)
        alpha:              (float) teleportation probability (default value 1)

    Returns:
        connected_disease_module:   list of genes containing the seed genes and the top-k
                                    ranked genes that form a connected component on the
                                    interactomne
    """

    d_entz_idx, d_idx_entz = create_mapping_index_entrezID(G)

    n_nodes = G.number_of_nodes()

    p0 = np.zeros(n_nodes)

    # select only the seed genes are on the PPI network
    seed_genes_on_PPI = [gene for gene in seed_genes if gene in d_entz_idx.keys()]

    # initialize (with optional scaling) of the visiting probability vector
    for gene in seed_genes_on_PPI:
        if scaling == 1:
            k = G.degree(gene)
            p0[d_entz_idx[gene]] = 1 * np.sqrt(k)
        else:
            p0[d_entz_idx[gene]] = 1.0

    # compute the colum-wise or symmetrical RW operator
    if symmetrical == 1:
        print("doing symmetrical")
        W = symmetric_rnd_walk_matrix(G, r=restart_parameter)
    else:
        print("doing column-wise")
        W = colwise_rnd_walk_matrix(G, r=restart_parameter, a=alpha)

    # apply the RW operator on the visiting probability vector (with optional scaling)
    if scaling == 1:
        Dinvsqrt = create_scaling_matrix(G)
        pinf = np.array(np.dot(Dinvsqrt, np.dot(W, p0)))
    else:
        pinf = np.dot(W, p0)

    # create dictionary of gene IDs and their corresponding visiting probability in sorted order
    d_gene_pvis_sorted = {}
    for p, x in sorted(zip(pinf, range(len(pinf))), reverse=True):
        d_gene_pvis_sorted[d_idx_entz[x]] = p / len(seed_genes_on_PPI)

    # obtain the ranking without seed genes
    rwr_rank_without_seed_genes = [g for g in d_gene_pvis_sorted.keys() if g not in seed_genes_on_PPI]

    # select the top ranked genes that lead to a connected component with the seed genes
    i = 0
    subgraph = nx.subgraph(G, seed_genes_on_PPI)
    connected_disease_module = [g for g in seed_genes_on_PPI]
    while not nx.is_connected(subgraph):
        connected_disease_module.append(rwr_rank_without_seed_genes[i])
        subgraph = nx.subgraph(G, connected_disease_module)
        i += 1

    with open(outfile_name, "w") as fout:
        fout.write("\t".join(["#rank", "RWR_node", "p_value"]))
        fout.write("\n")
        # fout.write('RWR_node' + '\t')
        rank = 0
        for g in connected_disease_module:
            rank += 1
            p = d_gene_pvis_sorted[g]
            fout.write("\t".join(map(str, ([rank, g, p]))))
            fout.write("\n")
            # fout.write(str(g) + '\t')

    return connected_disease_module


# =============================================================================

if __name__ == "__main__":
    input_list = sys.argv

    edgelist_file, seeds_file, scaling, sym, r, outfile_name = check_input_style(input_list)

    G, seed_genes = read_input(edgelist_file, seeds_file)

    connected_disease_module = rwr(G, seed_genes, scaling, sym, restart_parameter=r, alpha=1.0)
