process DRUGPREDICTIONS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(module)
    val idspace
    val algorithm
    val includeIndirectDrugs
    val includeNonApprovedDrugs
    val result_size

    output:
    tuple val(meta), path("${meta.id}.${algorithm}.drug_predictions.tsv")  , emit: drug_predictions
    tuple val(meta), path("${meta.id}.${algorithm}.csv"), emit: drugstone_download
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    drug_predictions.py --idspace "${idspace}" -p "${meta.id}" -a "${algorithm}" --includeIndirectDrugs ${includeIndirectDrugs} --includeNonApprovedDrugs ${includeNonApprovedDrugs} --result_size "${result_size}" "${module}" -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        drugstone: \$(pip show drugstone | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}
