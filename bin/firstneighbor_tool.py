#!/usr/bin/env python
from argparse import ArgumentParser
import graph_tool.all as gt


def main():
    """
    Extracts a subnetwork composed by the seeds and their first neighbors
    Execution example:
    python3 modulediscovery/bin/firstneighbor_tool.py -n modulediscovery-analysis/inputs/PPI.gt -s modulediscovery-analysis/inputs/seeds.txt -o modulediscovery-analysis/outputs/firstneighbor/subnetwork.gt
    """
    args = parse_user_arguments()
    run(args)


def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """
    description = 'First neighbor-based module identification'
    parser = ArgumentParser(description=description)
    parser.add_argument('-n', '--network_file', type=str, required=True, help='Path to file containing the network in graph-tool format')
    parser.add_argument('-s', '--seeds_file', type=str, required=True, help='Path to file containing the seeds')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to output file containing the resulting module in graph-tool format')
    args = parser.parse_args()
    return args


def run(args):
    """
    Runs the first neigbhor-based module identification
    """

    # Read the network
    g = gt.load_graph(str(args.network_file))

    # Read the seeds
    seeds = set()
    with open(args.seeds_file, "r") as seeds_fd:
        for line in seeds_fd:
            seeds.add(line.strip())

    # Iter over the nodes to get a dict mapping the IDs to their names
    node_name_to_id = {}
    for n in g.iter_vertices():
        node_name_to_id[g.vp["name"][n]] = n

    # Get the first neighbors of the seeds
    first_neighbors = set()
    for seed in seeds:
        if seed in node_name_to_id.keys():
            first_neighbors.add(node_name_to_id[seed])
            first_neighbors.update([n for n in g.iter_all_neighbors(node_name_to_id[seed])])

    # Create a property map to store a boolean values indicating whether the node is
    # part of the module (seed / neighbor) or not
    module_property = g.new_vertex_property("bool")

    # Mark the nodes based on if they are part of the module or not
    for node in g.vertices():
        if node in first_neighbors:
            module_property[node] = True
        else:
            module_property[node] = False

    # Extract the subgraph containing the module nodes
    subgraph = gt.GraphView(g, vfilt=module_property)

    # Save the subgraph in graph-tool format
    subgraph.save(args.output_file)

    return

if  __name__ == "__main__":
    main()
