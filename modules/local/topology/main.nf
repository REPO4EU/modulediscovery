process TOPOLOGY {
    tag "$meta.id"

    input:
    tuple val(meta), path(module)

    output:
    tuple val(meta), path("${meta.id}.topology.tsv")   , emit: topology
    path "versions.yml"                                , emit: versions

    script:
    """
    topology.py --module "$module" --id "${meta.id}" --out "${meta.id}.topology.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}


