#!/usr/bin/env python
import os
from argparse import ArgumentParser
import graph_tool.all as gt
import numpy as np


def main():
    """
    Filters a subnetwork based on the Subnetwork-Participation-Degree (SPD)
    Execution example:
    python3 modulediscovery/bin/spd_tool.py -s modulediscovery-analysis/outputs/diamond_seeds/firstneighbor/seeds_diamond.gt -n modulediscovery-analysis/inputs/PPI.gt -o modulediscovery-analysis/outputs/diamond_seeds/firstneighbor/subnetwork_filtered.gt -c 0.8
    """
    args = parse_user_arguments()
    run(args)


def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """
    description = "SPD-based module refinement"
    parser = ArgumentParser(description=description)
    parser.add_argument("-s", "--subnetwork_file", type=str, required=True, help="Path to file containing the subnetwork (disease module) in graph-tool format")
    parser.add_argument("-n", "--network_file", type=str, required=True, help="Path to file containing the network in graph-tool format")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="Path to output file containing the resulting module in graph-tool format")
    parser.add_argument("-c", "--cutoff", type=float, default=0.8, help="Fraction of the cumulative sum of SPD distribution values to be used as cut-off")
    args = parser.parse_args()
    return args


def run(args):
    """
    Runs a SPD-based module refinement and returns a filtered module
    """

    # Validate input files
    if not os.path.exists(args.subnetwork_file):
        raise FileNotFoundError(f"Subnetwork file not found: {args.subnetwork_file}")
    if not os.path.exists(args.network_file):
        raise FileNotFoundError(f"Network file not found: {args.network_file}")

    # Validate cutoff value
    if not (0 <= args.cutoff <= 1):
        raise ValueError("Cutoff must be a fraction between 0 and 1")

    # Read the subnetwork
    subnetwork = gt.load_graph(str(args.subnetwork_file))

    # Purge vertices without interactions
    subnetwork.purge_vertices()
    print(f"Subnetwork info: {subnetwork.num_vertices()} nodes and {subnetwork.num_edges()} edges")

    # Read the network
    full_interactome = gt.load_graph(str(args.network_file))
    print(f"Interactome info: {full_interactome.num_vertices()} nodes and {full_interactome.num_edges()} edges")

    # Check for the existence of the 'name' vertex property
    if 'name' not in subnetwork.vp:
        raise KeyError("Vertex property 'name' does not exist in the subnetwork graph")
    if 'name' not in full_interactome.vp:
        raise KeyError("Vertex property 'name' does not exist in the full interactome graph")

    # Create name to degree mappings for both networks
    name_to_degree_full = create_name_to_degree_map(full_interactome)
    name_to_degree_sub = create_name_to_degree_map(subnetwork)

    # Calculate SPD for each node in the subnetwork
    spd = subnetwork.new_vertex_property('float')
    for v in subnetwork.vertices():
        name = subnetwork.vp["name"][v]
        full_degree = name_to_degree_full.get(name, 0)
        sub_degree = name_to_degree_sub.get(name, 0)
        spd[v] = sub_degree / full_degree if full_degree > 0 else 0
    print(f"Max. SPD: {max(spd)}. Min. SPD: {min(spd)}")

    # Calculate the SPD cutoff
    sorted_spd = sorted(spd.a, reverse=True)
    cumulative_spd = np.cumsum(sorted_spd)
    cutoff_index = np.argmax(cumulative_spd >= cumulative_spd[-1] * args.cutoff)
    calculated_spd_cutoff = sorted_spd[cutoff_index]
    print(f"Calculated SPD cut-off: {calculated_spd_cutoff}")

    # Create a property map to store a boolean values indicating whether the node is
    # part of the pruned module or not
    module_property = subnetwork.new_vertex_property("bool")

    # Mark the nodes based on if they are part of the pruned module or not
    for node in subnetwork.vertices():
        if spd[node] >= calculated_spd_cutoff:
            module_property[node] = True
        else:
            module_property[node] = False

    # Extract the subgraph containing the pruned network nodes
    subnetwork_filtered = gt.GraphView(subnetwork, vfilt=module_property)
    print(f"Pruned subnetwork info: {subnetwork_filtered.num_vertices()} nodes and {subnetwork_filtered.num_edges()} edges")

    # Save the pruned network in graph-tool format
    subnetwork_filtered.save(args.output_file)

    return


#----------------------#
# Additional functions #
#----------------------#

def create_name_to_degree_map(graph):
    """
    Creates a dictionary mapping names to network degrees
    """
    degrees = graph.degree_property_map('total').a
    names = [graph.vp["name"][v] for v in graph.vertices()]
    return dict(zip(names, degrees))


if __name__ == "__main__":
    main()
