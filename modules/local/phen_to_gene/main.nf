process PHEN_TO_GENE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path (module)

    output:
    tuple val(meta), path("phenotype_to_gene.tsv"), emit: phen_gene
    path "versions.yml", emit: versions

    script:
    """
    prep_phen_to_gene.py --inpath ${module}

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
END_VERSIONS
    """
}
