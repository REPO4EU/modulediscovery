process FIRSTNEIGHBOR {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(seeds), path (network)

    output:
    tuple val(meta), path("${meta.id}.firstneighbor.gt"), emit: module
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    firstneighbor_tool.py -n $network -s $seeds -o "${meta.id}.firstneighbor.gt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
