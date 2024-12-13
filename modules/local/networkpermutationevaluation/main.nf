process NETWORKPERMUTATIONEVALUATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(module)
    path(permuted_modules)

    output:
    tuple val(meta), path("${meta.id}.network_permutation_multiqc_summary.tsv")     , emit: multiqc_summary
    tuple val(meta), path("${meta.id}.network_permutation_multiqc_jaccard.txt")     , emit: multiqc_jaccard
    path "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    network_permutation_evaluation.py \\
        --prefix ${meta.id} \\
        --module ${module} \\
        --permuted_modules ${permuted_modules} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
