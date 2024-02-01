process FIRSTNEIGHBOR {
    label 'process_single'

    conda 'conda-forge::graph-tool=2.58'
    container 'docker.io/tiagopeixoto/graph-tool:release-2.58'

    input:
    path seeds
    path network

    output:
    path "firstneighbor.gt", emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    firstneighbor_tool.py -n $network -s $seeds -o "firstneighbor.gt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
