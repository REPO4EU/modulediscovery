process GRAPHTOOLPARSER {
    label 'process_single'

    conda "conda-forge::graph-tool=2.58"
    container "docker.io/tiagopeixoto/graph-tool:release-2.58"

    input:
    path network
    val format

    output:
    path "*${format}*" , emit: network
    path "versions.yml", emit: versions

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
