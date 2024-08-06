process SAVEMODULES {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/graph-tool_networkx_pandas_pyvis_pruned:3b79bb1c134b59eb"

    input:
    tuple val(meta), path(module)

    output:
    tuple val(meta), path("${meta.id}.graphml")  , emit: graphml
    tuple val(meta), path("${meta.id}.nodes.tsv"), emit: nodes_tsv
    tuple val(meta), path("${meta.id}.edges.tsv"), emit: edges_tsv
    tuple val(meta), path("${meta.id}.pdf")      , emit: pdf
    tuple val(meta), path("${meta.id}.png")      , emit: png
    tuple val(meta), path("${meta.id}.svg")      , emit: svg
    tuple val(meta), path("${meta.id}.html")     , emit: html
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
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
        pyintergraph: \$(pip show pyintergraph | grep Version | awk '{print \$2}')
        pyvis: \$(pip show pyvis | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}
