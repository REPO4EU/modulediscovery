process DRUGPREDICTIONS {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/graph-tool_networkx_pandas_pyvis_pruned:3b79bb1c134b59eb"

    input:
    tuple val(meta), path(module)
    val idspace

    output:
    tuple val(meta), path("${meta.id}.drug_predictions.tsv")  , emit: drug_predictions
    tuple val(meta), path("${meta.id}.trustrank.csv"), emit: trustrank
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    drug_predictions.py --idspace "${idspace}" -p "${meta.id}" "${module}" -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        drugstone: \$(python -c "import drugstone; print(drugstone.__version__)")
    END_VERSIONS
    """
}
