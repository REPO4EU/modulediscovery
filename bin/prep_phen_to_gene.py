#!/usr/bin/env python3

import os
import glob
import pandas
import argparse


def map_n_restructure(inpath, outpath, id_mapping_file=None):
    list_files = glob.glob(os.path.join(inpath, "modules/tsv_nodes/*.nodes.tsv"))
    dict_genes = {}

    if id_mapping_file:
        with open(id_mapping_file) as mapping:
            mapping = {k: gene for k, gene in (l.strip().split("\t") for l in mapping)}
    for f in list_files:
        base = os.path.basename(f)
        name = base.replace(".nodes.tsv", "")  # name to use for .gmt file
        # Read content.
        file_content = pandas.read_csv(f, sep="\t")
        genes = list(file_content["name"])

        if id_mapping_file:  # Perform mapping only if mapping file is provided
            mapped_genes = [mapping.get(str(gene), gene) for gene in genes]
        else:
            mapped_genes = genes

        dict_genes[name] = mapped_genes
    with open(os.path.join(outpath + "phenotype_to_gene.tsv"), "w") as fout:
        fout.write("phenotype\tgene\n")  # write header
        for key, list_values in dict_genes.items():
            for gene_val in list_values:
                line = f"{key}\t{gene_val}\n"
                fout.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map enrez IDs to gene symbols and prepare phenotype to gene file."
    )
    parser.add_argument(
        "--inpath",
        type=str,
        default="/Users/htoukabri/myFolder/testdisc/testout/",
        help="Path to modules folder.",
    )
    parser.add_argument(
        "--outpath", type=str, default="./", help="Path to output folder."
    )
    parser.add_argument(
        "--id_mapping_file",
        type=str,
        default=None,
        help="Path to the ID-gene symbols mapping file.",
    )

    args = parser.parse_args()
    map_n_restructure(args.inpath, args.outpath, args.id_mapping_file)
