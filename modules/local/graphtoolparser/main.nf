process GRAPHTOOLPARSER {
    label 'process_single'

    input:
    path network
    val format

    output:
    path "*${format}*"          , emit: network
    path "input_network_mqc.tsv", emit: multiqc, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    graph_tool_parser.py $network -f $format -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
