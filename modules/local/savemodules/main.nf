process SAVEMODULES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(module)

    output:
    tuple val(meta), path("${meta.id}.graphml")  , emit: graphml
    tuple val(meta), path("${meta.id}.nodes.tsv"), emit: nodes_tsv
    tuple val(meta), path("${meta.id}.edges.tsv"), emit: edges_tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    save_modules.py -m "${module}" -p "${meta.id}" -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
