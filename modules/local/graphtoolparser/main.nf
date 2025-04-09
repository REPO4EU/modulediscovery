process GRAPHTOOLPARSER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), (path(network), stageAs: 'input/*')
    val format

    output:
    tuple val(meta), path("*${format}*")               , emit: network
    tuple val(meta), path("input_network_multiqc.tsv") , emit: multiqc, optional: true
    path "versions.yml"                                , emit: versions

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
