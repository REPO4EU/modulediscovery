process PROXIMITY {
    label 'process_single'

    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    path network
    path shortest_paths
    path drug_to_target
    path phen_to_gene

    output:
    path("*.tsv"), emit: proxout
    path "versions.yml", emit: versions

    script:
    def shortest_paths = shortest_paths.name != 'NO_FILE' ? "$shortest_paths" : 'None'
    def drug_to_target = drug_to_target.name != 'NO_FILE' ? "$drug_to_target" : 'None'
    """
    # Create a config file.
    cat <<EOT > proximity_config.txt
    [PROXIMITY]
    drug_to_target = ${drug_to_target}
    drug_column = drugbankId
    target_column = targetDomainId
    phenotype_to_gene = ${phen_to_gene}
    phenotype_column = phenotype
    gene_column = gene
    network_file = ${network}
    shortest_paths = ${shortest_paths}
    id_mapping_file = None
    output_file = proximity.tsv
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
