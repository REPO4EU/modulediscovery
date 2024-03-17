#! /usr/bin/env python

import argparse
import os
import pandas as pd
from biodigest.single_validation import single_validation, save_results
from biodigest.evaluation.d_utils.plotting_utils import create_plots,create_extended_plots

def run_analysis (target_set, target_type, network, network_type, output_directory):
    # ==== define required input ====
    tar_set = pd.read_csv(target_set, header=None, sep="\t", dtype=str)[0]
    tar_id_type = target_type
    mode = "subnetwork"

    # ==== define input for network integration ====
    network_data = {"network_file":network, "prop_name":"name", "id_type": network_type}

    # ==== define optional input influencing results ====
    distance_measure="jaccard"
    background_model="network" # for subnetwork mode the only option so it is set fixed
    runs = 1000 # how many random runs for empirical p-value estimation
    perc = 100 # how many % of the original input should be perturbated for the background model
    enriched=False

    # ==== define optional input influencing saving of results ====
    out_dir = output_directory
    verbose=True # printing additional information during the run
    prefix="subnetwork"

    results = single_validation(tar=tar_set, tar_id=tar_id_type, mode=mode, runs=runs, background_model=background_model, verbose=verbose, enriched=enriched, distance=distance_measure, network_data=network_data)

    pd.DataFrame(results["p_values"]['values'])

    pd.DataFrame(results["input_values"]['values'])

    save_results(results=results, prefix=prefix, out_dir=out_dir)

    create_plots(results=results, mode=mode, tar=tar_set, tar_id=tar_id_type, out_dir=out_dir, prefix=prefix, file_type='png')

    create_extended_plots(results=results, mode=mode, tar=tar_set, out_dir=out_dir, prefix=prefix, file_type='png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Execute digest")

    parser.add_argument("--target_file", required=True, help="Input target file path")
    parser.add_argument("--target_type", required=True, help="Input target id")
    parser.add_argument("--network", required=True, help="Input network")
    parser.add_argument("--network_type", required=True, help="Input network id")

    args = parser.parse_args()
    target_file = args.target_file
    target_type = args.target_type
    network = args.network
    network_type = args.network_type   
    output_directory = os.path.splitext(os.path.basename(target_file))[0]
    os.makedirs(output_directory, exist_ok=True)

    run_analysis(target_file, target_type, network, network_type, output_directory)