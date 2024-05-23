#!/usr/bin/env python
import os
from argparse import ArgumentParser
import graph_tool.all as gt
import numpy as np
import scipy


def main():
    """
    Filters a subnetwork based on the Subnetwork-Participation-Degree (SPD)
    Execution examples:
    python3 modulediscovery/bin/spd_filter_tool.py -s modulediscovery-analysis/outputs/thyroid_cancer_intogen/firstneighbor/firstneighbor.gt -n modulediscovery-analysis/outputs/thyroid_cancer_intogen/graphtoolparser/nedrex_ppi_genename_20240205_nedrex.gt -o modulediscovery-analysis/outputs/thyroid_cancer_intogen/firstneighbor/firstneighbor_spd.gt -t fraction -c 0.95
    python3 modulediscovery/bin/spd_filter_tool.py -s modulediscovery-analysis/outputs/thyroid_cancer_tcga_threshold_200/firstneighbor/firstneighbor.gt -n modulediscovery-analysis/outputs/thyroid_cancer_tcga_threshold_200/graphtoolparser/nedrex_ppi_genename_20240205_nedrex.gt -o modulediscovery-analysis/outputs/thyroid_cancer_tcga_threshold_200/firstneighbor/firstneighbor_spd.gt -t zscore_fraction -c 0.95
    """
    args = parse_user_arguments()
    run(args)


def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """
    description = "SPD-based module refinement"
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
    parser.add_argument(
        "-t",
        "--type_cutoff",
        type=str,
        default="fraction",
        help="Type of cut-off. Options: fraction, zscore_fraction, spd, spd_mean_distribution_elbow",
    )
    parser.add_argument(
        "-c", "--cutoff", type=float, default=0.95, help="Cut-off threshold"
    )
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

    # Validate cutoff type
    valid_type_cutoffs = [
        "fraction",
        "zscore_fraction",
        "spd",
        "spd_mean_distribution_elbow",
    ]
    if args.type_cutoff not in valid_type_cutoffs:
        raise ValueError(
            "Invalid input. Please enter one of these: fraction, zscore_fraction, spd, spd_mean_distribution_elbow"
        )
    print(f"SPD parameters info: type {args.type_cutoff} and value {args.cutoff}")

    # Validate cutoff value
    if args.type_cutoff != "spd_mean_distribution_elbow" and not (
        0 <= args.cutoff <= 1
    ):
        raise ValueError("Cutoff must be a fraction between 0 and 1")

    # Read the subnetwork
    subnetwork = gt.load_graph(str(args.subnetwork_file))

    # Purge vertices without interactions
    subnetwork.purge_vertices()
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
    spd_sorted = sorted(spd.a, reverse=True)
    print(f"Max. SPD: {max(spd)}. Min. SPD: {min(spd)}")

    # Calculate SPD cut-off based on the fraction of nodes in the SPD distribution
    if args.type_cutoff == "fraction":
        rank_cutoff = round(len(spd_sorted) * (1 - args.cutoff))
        # If there are multiple ranks with the same SPD, we consider all of them, even if we end up
        # having a fraction of nodes larger than the initially considered
        spd_cutoff = spd_sorted[rank_cutoff - 1]

    # Calculate SPD cut-off based on the fraction of nodes in the z-score SPD distribution
    elif args.type_cutoff == "zscore_fraction":
        spd_mean = np.mean(list(spd.a))
        spd_std = np.std(list(spd.a))
        spd_zscore = [(spd_val - spd_mean) / spd_std for spd_val in spd.a]
        spd_zscore.sort()
        z_cutoff = scipy.stats.norm.ppf(
            args.cutoff, loc=np.mean(spd_zscore), scale=np.std(spd_zscore)
        )
        spd_cutoff = z_cutoff * spd_std + spd_mean
        if spd_cutoff > 1:
            print(
                f"Calculated SPD cut-off higher than 1: {round(spd_cutoff, 3)}. We will use as cut-off a SPD cut-off of 1 instead."
            )
            spd_cutoff = 1

    # Use SPD cut-off
    elif args.type_cutoff == "spd":
        spd_cutoff = args.cutoff

    # Calculate SPD cut-off based on the elbow of the mean SPD distribution
    elif args.type_cutoff == "spd_mean_distribution_elbow":
        (mean_spd, _) = calculate_mean_spd_distribution(
            subnetwork, spd, name_to_degree_full
        )
        spd_rank = list(range(1, len(mean_spd) + 1))
        rank_cutoff = find_elbow_point(x=spd_rank, y=mean_spd, plot=False)
        print(f"Elbow rank cut-off: {rank_cutoff}")
        spd_cutoff = spd_sorted[rank_cutoff - 1]

    print(f"SPD cut-off: {spd_cutoff}")

    # Create a property map to store a boolean values indicating whether the node is
    # part of the pruned module or not
    module_property = subnetwork.new_vertex_property("bool")

    # Mark the nodes based on if they are part of the pruned module or not
    for node in subnetwork.vertices():
        if spd[node] >= spd_cutoff:
            module_property[node] = True
        else:
            module_property[node] = False

    # Extract the subgraph containing the pruned network nodes
    subnetwork_filtered = gt.GraphView(subnetwork, vfilt=module_property)
    print(
        f"Pruned subnetwork info: {subnetwork_filtered.num_vertices()} nodes and {subnetwork_filtered.num_edges()} edges"
    )

    # Save the pruned network in graph-tool format
    subnetwork_filtered.save(args.output_file)

    return


# ----------------------#
# Additional functions #
# ----------------------#


def create_name_to_degree_map(graph):
    """
    Creates a dictionary mapping names to network degrees
    """
    degrees = graph.degree_property_map("total").a
    names = [graph.vp["name"][v] for v in graph.vertices()]
    return dict(zip(names, degrees))


def calculate_spd_subnetwork(subnetwork, name_to_degree_sub, name_to_degree_full):
    """
    Calculates the spd of all the nodes in a subnetwork.
    Returns the spd in form of graph_tool vertex property and the subnetwork containing
    the spd as vertex property.
    """
    # spd = subnetwork.new_vertex_property('float')
    subnetwork.vp["spd"] = subnetwork.new_vertex_property("float")
    names = subnetwork.vp["name"]
    full_degrees = np.array([name_to_degree_full.get(name, 0) for name in names])
    sub_degrees = np.array([name_to_degree_sub.get(name, 0) for name in names])

    # Avoid division by zero and calculate SPD
    with np.errstate(divide="ignore", invalid="ignore"):
        spd_values = np.true_divide(sub_degrees, full_degrees)
        spd_values[full_degrees == 0] = 0  # Set SPD to 0 where full_degree is 0

    # Check for SPD values greater than 1
    if np.any(spd_values > 1):
        raise Exception("ERROR: Node with SPD higher than 1 detected.")

    # Assign the values back to the graph-tool property map
    subnetwork.vp["spd"].get_array()[:] = spd_values

    return subnetwork.vp["spd"], subnetwork


def calculate_mean_spd_distribution(
    subnetwork, vertex_property_to_sort, name_to_degree_full
):
    """
    Starting from the nodes with higher SPD or score (given by vertex_property_to_sort),
    iteratively includes nodes to a subnetwork and calculates the mean SPD for each
    subnetwork.
    """

    # Sort nodes based on SPD from high to low
    sorted_nodes = sorted(
        subnetwork.vertices(), key=lambda v: vertex_property_to_sort[v], reverse=True
    )

    # Initialize variables for the subsubnetworks
    subsubnetwork_property = subnetwork.new_vertex_property("bool")
    mean_spd_list = []
    median_spd_list = []

    # Add new node to subnetwork, recalculate SPD, and calculate mean SPD
    for v in sorted_nodes:

        # Add node to subsubnetwork
        subsubnetwork_property[v] = True

        # Extract the subgraph containing the subsubnetwork nodes
        subsubnetwork = gt.GraphView(subnetwork, vfilt=subsubnetwork_property)
        subsubnetwork.purge_vertices()

        # Create name to degree mappings for the subsubnetwork
        name_to_degree_subsub = create_name_to_degree_map(subsubnetwork)

        # Recalculate SPD for each node in the subsubnetwork
        subsub_spd, subsubnetwork = calculate_spd_subnetwork(
            subsubnetwork, name_to_degree_subsub, name_to_degree_full
        )

        # Calculate mean SPD so far
        mean_spd_list.append(np.mean(list(subsub_spd)))
        median_spd_list.append(np.median(list(subsub_spd)))
        # print(f"Node ID {v} and name {subsubnetwork.vp["name"][v]}. SPD: {vertex_property_to_sort[v]}. Num nodes in subnetwork: {subsubnetwork.num_vertices()}. Mean SPD: {np.mean(list(subsub_spd))}. Median SPD: {np.median(list(subsub_spd))}")

    return mean_spd_list, median_spd_list


def find_elbow_point(x, y, plot=False):
    """
    Find the elbow point in a curve defined by x and y.

    Parameters:
    - x: A list or array of x values.
    - y: A list or array of y values.
    - plot: Boolean, whether to plot the curve and the elbow point.

    Returns:
    - elbow_index: The index of the elbow point in the x and y arrays.
    """
    # Convert lists to numpy arrays if they aren't already
    x = np.array(x)
    y = np.array(y)

    # Define the start and end points of the curve
    start_point = np.array([x[0], y[0]])
    end_point = np.array([x[-1], y[-1]])

    # Calculate a line vector from start to end
    line_vec = end_point - start_point

    # Calculate vectors from the start point to each point in the curve
    point_vecs = np.column_stack((x, y)) - start_point

    # Normalize the line vector to create a unit vector representing the direction of
    # the line
    line_unit_vec = line_vec / np.linalg.norm(line_vec)

    # Project each point vector onto the normalized line vector. This projection gives
    # the distance of each point from the line connecting the start and end points
    proj_lengths = np.dot(point_vecs, line_unit_vec)

    # Calculate the actual points on the line
    proj_points = np.outer(proj_lengths, line_unit_vec)

    # Calculate the distance from the points to the line
    distances = np.linalg.norm(point_vecs - proj_points, axis=1)

    # Find the index of the point with the maximum distance
    elbow_index = np.argmax(distances)

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(x, y, label="Curve")
        ax.scatter(x[elbow_index], y[elbow_index], color="red", label="Elbow Point")
        ax.legend()
        plt.show()
        return elbow_index, fig
    else:
        return elbow_index


if __name__ == "__main__":
    main()
