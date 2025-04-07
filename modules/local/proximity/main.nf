process PROXIMITY {
    tag "$meta.id"
    label 'process_single'

    input:
    path network
    path shortest_paths
    path drug_to_target
    tuple val(meta), path (module)

    output:
    path("${meta.id}.proximity.tsv"), emit: proxout
    path "versions.yml", emit: versions

    script:
    """
    # Create a config file.
    cat <<EOT > proximity_config.txt
    [PROXIMITY]
    drug_to_target = ${drug_to_target}
    drug_column = drugbankId
    target_column = targetDomainId
    phenotype_to_gene = ${module}
    phenotype_column = None
    prefix = ${meta.id}
    gene_column = name
    network_file = ${network}
    shortest_paths = ${shortest_paths}
    id_mapping_file = None
    output_file = ${meta.id}.proximity.tsv
    EOT

    # Run proximity.
    proximity.py proximity_config.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
    END_VERSIONS
    """
}
