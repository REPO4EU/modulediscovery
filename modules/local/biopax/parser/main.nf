process BIOPAX_PARSER {
    label 'process_single'

//     conda "conda-forge::graph-tool=2.58"
    container "quirinmanz/gt2biopax:latest" //TODO: should be fixed tag

    input:
    path network
    val format

    output:
    path "*.owl" , emit: biopax
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gt2biopax.py $network -i $format -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
        pybiopax: \$(python -c "import pybiopax; print(pybiopax.__version__)")
    END_VERSIONS
    """
}
