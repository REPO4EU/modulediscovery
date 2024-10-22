process PHEN_TO_GENE {
    label 'process_single'

    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    path modules

    output:
    path("*.tsv"), emit: phen_gene
    path "versions.yml", emit: versions

    script:
    """
    prep_phen_to_gene.py --inpath ${modules}

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
END_VERSIONS
    """
}
