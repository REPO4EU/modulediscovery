process BIOPAX_PARSER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(network)
    val idspace
    val add_variants

    output:
    path "*.owl" , emit: biopax
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gt2biopax.py $network -i $idspace -l DEBUG ${add_variants ? '-v' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
        pybiopax: \$(python -c "import pybiopax; print(pybiopax.__version__)")
        nedrex: \$(python -c "import nedrex; print(nedrex.__version__)")
    END_VERSIONS
    """
}
