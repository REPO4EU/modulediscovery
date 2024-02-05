process MODULEPARSER {
    label 'process_single'

    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    path network
    val format
    path module

    output:
    path "${format}.gt" , emit: network
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    graph_tool_parser.py $network -f $format -l DEBUG -m $module

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
