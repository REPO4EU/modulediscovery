#!/usr/bin/env python
import os
from argparse import ArgumentParser
import graph_tool.all as gt
import numpy as np


def main():
    """
    Annotates the Subnetwork-Participation-Degree (SPD) and component ID in a subnetwork
    Execution example:
    python3 modulediscovery/bin/network_annotation.py -s modulediscovery-analysis/outputs/thyroid_cancer_intogen/firstneighbor/firstneighbor.gt -n modulediscovery-analysis/outputs/thyroid_cancer_intogen/graphtoolparser/nedrex_ppi_genename_20240205_nedrex.gt -o modulediscovery-analysis/outputs/thyroid_cancer_intogen/firstneighbor/firstneighbor_spd.gt
    """
    args = parse_user_arguments()
    run(args)


def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """
    description = "Network annotation of SPD and component ID"
    parser = ArgumentParser(description=description)
    parser.add_argument(
        "-s",
        "--subnetwork_file",
        type=str,
        required=True,
        help="Path to file containing the subnetwork (disease module) in graph-tool format",
    )
    parser.add_argument(
        "-n",
        "--network_file",
        type=str,
        required=True,
        help="Path to file containing the network in graph-tool format",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="Path to output file containing the resulting module in graph-tool format",
    )
    args = parser.parse_args()
    return args


def run(args):
    """
    Calculates the SPD of the nodes in a module
    """

    # Validate input files
    if not os.path.exists(args.subnetwork_file):
        raise FileNotFoundError(f"Subnetwork file not found: {args.subnetwork_file}")
    if not os.path.exists(args.network_file):
        raise FileNotFoundError(f"Network file not found: {args.network_file}")

    # Read the subnetwork
    subnetwork = gt.load_graph(str(args.subnetwork_file))

    # Purge vertices without interactions
    print(
        f"Subnetwork info: {subnetwork.num_vertices()} nodes and {subnetwork.num_edges()} edges"
    )

    # Read the network
    full_interactome = gt.load_graph(str(args.network_file))
    print(
        f"Interactome info: {full_interactome.num_vertices()} nodes and {full_interactome.num_edges()} edges"
    )

    # Check for the existence of the 'name' vertex property in the graph_tool networks
    if "name" not in subnetwork.vp:
        raise KeyError("Vertex property 'name' does not exist in the subnetwork graph")
    if "name" not in full_interactome.vp:
        raise KeyError(
            "Vertex property 'name' does not exist in the full interactome graph"
        )

    # Create name to degree mappings for both networks
    name_to_degree_full = create_name_to_degree_map(full_interactome)
    name_to_degree_sub = create_name_to_degree_map(subnetwork)

    # Calculate SPD for each node in the subnetwork
    spd, subnetwork = calculate_spd_subnetwork(
        subnetwork, name_to_degree_sub, name_to_degree_full
    )

    # Assign component ID to each component of the subnetwork
    component_id = assign_component_ids(subnetwork)

    # Save the network containing the annotations in graph-tool format
    subnetwork.save(args.output_file)

    return


# ----------------------#
# Additional functions #
# ----------------------#


def create_name_to_degree_map(graph):
    """
    Creates a dictionary mapping names to network degrees
    """
    name_to_degree = {graph.vp["name"][v]: v.out_degree() for v in graph.vertices()}
    return name_to_degree


def calculate_spd_subnetwork(subnetwork, name_to_degree_sub, name_to_degree_full):
    """
    Calculates the spd of all the nodes in a subnetwork.
    Returns the spd in form of graph_tool vertex property and the subnetwork containing
    the spd as vertex property.
    """
    subnetwork.vp["spd"] = subnetwork.new_vertex_property("float")
    for v in subnetwork.vertices():
        name = subnetwork.vp["name"][v]
        full_degree = name_to_degree_full.get(name, 0)
        sub_degree = name_to_degree_sub.get(name, 0)
        subnetwork.vp["spd"][v] = sub_degree / full_degree if full_degree > 0 else 0

    # Check for SPD values greater than 1
    if np.any(subnetwork.vp["spd"].a > 1):
        raise Exception("ERROR: Node with SPD higher than 1 detected.")

    return subnetwork.vp["spd"], subnetwork


def assign_component_ids(graph):
    graph.vp["component_id"] = graph.new_vertex_property("int")
    component_property_map, component_histogram = gt.label_components(graph)

    for v in graph.vertices():
        graph.vp["component_id"][v] = component_property_map[v]

    return graph.vp["component_id"]


if __name__ == "__main__":
    main()
