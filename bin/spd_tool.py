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
    spd = calculate_spd_subnetwork(subnetwork, name_to_degree_sub, name_to_degree_full)
    print(f"Max. SPD: {max(spd)}. Min. SPD: {min(spd)}")

    ## Calculate the SPD average for each subnetwork, starting from the top-ranked nodes
    # Sort nodes based on SPD from high to low
    sorted_nodes = sorted(subnetwork.vertices(), key=lambda v: spd[v], reverse=True)

    # Initialize variables for the subsubnetworks
    subsubnetwork_property = subnetwork.new_vertex_property("bool")
    mean_spd_list = []
    median_spd_list = []

    # Iterate, add to subnetwork, recalculate SPD, and calculate mean SPD
    for v in sorted_nodes:

        # Add node to subsubnetwork
        subsubnetwork_property[v] = True

        # Extract the subgraph containing the subsubnetwork nodes
        subsubnetwork = gt.GraphView(subnetwork, vfilt=subsubnetwork_property)
        subsubnetwork.purge_vertices()

        # Create name to degree mappings for the subsubnetwork
        name_to_degree_subsub = create_name_to_degree_map(subsubnetwork)

        # Recalculate SPD for each node in the subsubnetwork
        subsub_spd = calculate_spd_subnetwork(subsubnetwork, name_to_degree_subsub, name_to_degree_full)

        # Calculate mean SPD so far
        mean_spd_list.append(np.mean(list(subsub_spd)))
        median_spd_list.append(np.median(list(subsub_spd)))
        #print(f"Node ID {v} and name {subsubnetwork.vp["name"][v]}. SPD: {spd[v]}. Num nodes in subnetwork: {subsubnetwork.num_vertices()}. Mean SPD: {np.mean(list(subsub_spd))}. Median SPD: {np.median(list(subsub_spd))}")

    print(mean_spd_list)
    #print(median_spd_list)
    print(np.mean(mean_spd_list))

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
    names = graph.vp["name"]
    return dict(zip(names, degrees))

def calculate_spd_subnetwork(subnetwork, name_to_degree_sub, name_to_degree_full):
    """
    Calculates the spd of all the nodes in a subnetwork.
    Returns the spd in form of graph_tool vertex property.
    """

    spd = subnetwork.new_vertex_property('float')
    names = subnetwork.vp["name"]
    full_degrees = np.array([name_to_degree_full.get(name, 0) for name in names])
    sub_degrees = np.array([name_to_degree_sub.get(name, 0) for name in names])

    # Avoid division by zero and calculate SPD
    with np.errstate(divide='ignore', invalid='ignore'):
        spd_values = np.true_divide(sub_degrees, full_degrees)
        spd_values[full_degrees == 0] = 0  # Set SPD to 0 where full_degree is 0

    # Check for SPD values greater than 1
    if np.any(spd_values > 1):
        raise Exception("ERROR: Node with SPD higher than 1 detected.")

    # Assign the values back to the graph-tool property map
    spd.get_array()[:] = spd_values

    return spd

if __name__ == "__main__":
    main()
