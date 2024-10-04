process PROXIMITY {
    tag "$meta.id"
    label 'process_single'

    input:
    path network
    tuple val(meta), path(modules)
    val shortest_paths
    val drug_to_target

    output:
    path("*.tsv"), emit: proxout
    path("*.pickle"), optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Reformat the disease module output.
    prep_phen_to_gene.py --inpath ${launchDir}/${params.outdir}/ 

    # Create a temporary configuration file.
    config_file=\$(mktemp)

cat <<EOT > \$config_file
    [PROXIMITY]
    drug_to_target = $drug_to_target
    drug_column = drugbankId
    target_column = targetDomainId
    phenotype_to_gene = phenotype_to_gene.tsv
    phenotype_column = phenotype
    gene_column = gene
    network_file = $network 
    shortest_paths = $shortest_paths
    id_mapping_file = None
    output_file = proximity.tsv
EOT

    # Run proximity.
    proximity.py \$config_file 

    # Delete the configfile once proximity script is done.
    rm \$config_file

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
END_VERSIONS
    """
}

