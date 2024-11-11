#!/usr/bin/env python

"""Randomize an input network file (degree-preserving randomization)."""

import argparse
import logging
import sys
from pathlib import Path
import networkx as nx
import math, random
import graph_tool.all as gt
import util
import pyintergraph
import copy


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
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def connected_double_edge_swap(
    G, is_weighted=False, nswap=5, _window_threshold=3, seed=random
):
    """Attempts the specified number of double-edge swaps in the graph `G`.

        A double-edge swap removes two randomly chosen edges `(u, v)` and `(x,
        y)` and creates the new edges `(u, x)` and `(v, y)`:

         u--v            u  v
                        becomes  |  |
         x--y            x  y

        If either `(u, x)` or `(v, y)` already exist, then no swap is performed
        so the actual number of swapped edges is always *at most* `nswap`.

        Parameters
        ----------
        G : networkx graph
           An undirected graph

    is_weighted : boolean (default = False)
       Consider a weighted network or not

        nswap : integer (optional, default=5)
           Number of double-edge swaps to perform

        _window_threshold : integer (default = 3)
           The window size below which connectedness of the graph will be checked
           after each swap.

           The "window" in this function is a dynamically updated integer that
           represents the number of swap attempts to make before checking if the
           graph remains connected. It is an optimization used to decrease the
           running time of the algorithm in exchange for increased complexity of
           implementation.

           If the window size is below this threshold, then the algorithm checks
           after each swap if the graph remains connected by checking if there is a
           path joining the two nodes whose edge was just removed. If the window
           size is above this threshold, then the algorithm performs do all the
           swaps in the window and only then check if the graph is still connected.

        seed : integer, random (default), or None
                Indicator of random number generation state.
                See :ref:`Randomness<randomness>`.

        Returns
        -------
    networkx graph : the randomized network

        int :  The number of successful swaps

        Raises
        ------

        NetworkXError

           If the input graph is not connected, or if the graph has fewer than four
           nodes.

        Notes
        -----

        The initial graph `G` must be connected, and the resulting graph is
        connected. The graph `G` is modified in place.

        References
        ----------
        .. [1] C. Gkantsidis and M. Mihail and E. Zegura,
                   The Markov chain simulation method for generating connected
                   power law random graphs, 2003.
                   https://faculty.cc.gatech.edu/~ewz/papers/markovgen.pdf
    """
    if not nx.is_connected(G):
        raise nx.NetworkXError("Graph not connected")
    if len(G) < 4:
        raise nx.NetworkXError("Graph has less than four nodes.")
    n = 0
    swapcount = 0
    # Label key for nodes
    dk = list(n for n, d in G.degree())
    cdf = nx.utils.cumulative_distribution(list(d for n, d in G.degree()))
    discrete_sequence = nx.utils.discrete_sequence
    window = 1

    while n < nswap:
        wcount = 0
        swapped = []
        # If the window is small, we just check each time whether the graph is
        # connected by checking if the nodes that were just separated are still
        # connected.
        if window < _window_threshold:
            # This Boolean keeps track of whether there was a failure or not.
            fail = False
            while wcount < window and n < nswap:
                # Pick two random edges without creating the edge list. Choose
                # source nodes from the discrete degree distribution.
                (ui, xi) = discrete_sequence(2, cdistribution=cdf, seed=seed)
                # If the source nodes are the same, skip this pair.
                while ui == xi:
                    (ui, xi) = nx.utils.discrete_sequence(2, cdistribution=cdf)
                # Convert an index to a node label.
                u = dk[ui]
                x = dk[xi]
                # Choose targets uniformly from neighbors.
                v = seed.choice(list(G.neighbors(u)))
                y = seed.choice(list(G.neighbors(x)))

                # If the target nodes are the same, skip this pair.
                if v == y:
                    continue

                if x not in G[u] and y not in G[v]:
                    if is_weighted:
                        G.add_edge(u, x, weight=G[u][v]["weight"])
                        G.add_edge(v, y, weight=G[x][y]["weight"])
                    else:
                        G.add_edge(u, x)
                        G.add_edge(v, y)
                    G.remove_edge(u, v)
                    G.remove_edge(x, y)

                    swapped.append((u, v, x, y))
                    swapcount += 1
                else:
                    fail = True
                n += 1

                # If G remains connected...
                if nx.has_path(G, u, v):
                    wcount += 1

                # Otherwise, undo the changes.
                else:
                    if is_weighted:
                        G.add_edge(u, v, weight=G[u][x]["weight"])
                        G.add_edge(x, y, weight=G[v][y]["weight"])
                    else:
                        G.add_edge(u, v)
                        G.add_edge(x, y)
                    G.remove_edge(u, x)
                    G.remove_edge(v, y)
                    swapcount -= 1
                    fail = True
            # If one of the swaps failed, reduce the window size.
            if fail:
                window = int(math.ceil(window / 2))
            else:
                window += 1
        # If the window is large, then there is a good chance that a bunch of
        # swaps will work. It's quicker to do all those swaps first and then
        # check if the graph remains connected.
        else:
            while wcount < window and n < nswap:
                # Pick two random edges without creating the edge list. Choose
                # source nodes from the discrete degree distribution.
                (ui, xi) = nx.utils.discrete_sequence(2, cdistribution=cdf)
                # If the source nodes are the same, skip this pair.
                while ui == xi:
                    (ui, xi) = nx.utils.discrete_sequence(2, cdistribution=cdf)

                # Convert an index to a node label.
                u = dk[ui]
                x = dk[xi]
                # Choose targets uniformly from neighbors.
                v = seed.choice(list(G.neighbors(u)))
                y = seed.choice(list(G.neighbors(x)))

                # If the target nodes are the same, skip this pair.
                if v == y:
                    continue

                if x not in G[u] and y not in G[v]:
                    if is_weighted:
                        G.add_edge(u, x, weight=G[u][v]["weight"])
                        G.add_edge(v, y, weight=G[x][y]["weight"])
                    else:
                        G.add_edge(u, x)
                        G.add_edge(v, y)
                    G.remove_edge(u, v)
                    G.remove_edge(x, y)
                    swapped.append((u, v, x, y))
                    swapcount += 1
                n += 1
                wcount += 1

            # If G remains connected and there is no de-swaps, increase the window size
            if nx.is_connected(G):
                window += 1

            # Otherwise, undo the changes from the previous window and decrease
            # the window size.
            else:
                removed_edges = []
                while swapped:
                    (u, v, x, y) = swapped.pop()
                    if is_weighted:
                        G.add_edge(u, v, weight=G[u][x]["weight"])
                        G.add_edge(x, y, weight=G[v][y]["weight"])
                    else:
                        G.add_edge(u, v)
                        G.add_edge(x, y)

                    try:
                        G.remove_edge(u, x)
                        removed_edges.append(sorted((u, x)))
                    except:
                        print(("PROBLEM"))
                        print("failed: ", u, x)
                        print(x in G[u], u in G[x])
                        print(G[u])
                        print("already removed edges:", removed_edges)

                    try:
                        G.remove_edge(v, y)
                        removed_edges.append(sorted((v, y)))
                    except:
                        print(("PROBLEM"))
                        print("already removed: ", u, x, "failed: ", v, y)
                        print(v in G[y], y in G[v])
                        print(G[v])
                        print("already removed edges:", removed_edges)

                    swapcount -= 1
                window = int(math.ceil(window / 2))

    return [G, swapcount]


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.network.is_file():
        logger.error(f"The given input file {args.network} was not found!")
        sys.exit(2)
    logger.debug(f"{args=}")

    stem = Path(args.network).stem
    extension = Path(args.network).suffix

    g = util.load_graph(str(args.network))
    nx_graph = pyintergraph.gt2nx(g, labelname="name")

    weighted = nx.is_weighted(nx_graph)

    if not nx.is_connected(nx_graph):

        logging.warning(
            "Provided network is not connected! \
                        Will output randomized network from largest connected component."
        )
        gcc = sorted(nx.connected_components(nx_graph), key=len, reverse=True)
        nx_connected_graph = nx_graph.subgraph(gcc[0]).copy()
        print(nx.is_frozen(nx_connected_graph))

    else:
        nx_connected_graph = copy.deepcopy(nx_graph)

    for i in range(1000):
        nx_G = connected_double_edge_swap(nx_connected_graph, is_weighted=weighted)[0]
        gt_G = pyintergraph.nx2gt(nx_G, labelname="name")
        output_file = f"{stem}.perm_{i}{extension}"
        gt_G.save(output_file)


if __name__ == "__main__":
    sys.exit(main())
