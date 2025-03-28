
process SEEDPERMUTATIONEVALUATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(module)
    path(seeds)
    path(permuted_modules)
    path(permuted_seeds)
    path(network)

    output:
    tuple val(meta), path("${meta.id}.seed_permutation_evaluation_summary.tsv")
    tuple val(meta), path("${meta.id}.seed_permutation_evaluation_detailed.tsv")
    tuple val(meta), path("${meta.id}.seed_permutation_multiqc_summary.tsv")     , emit: multiqc_summary
    tuple val(meta), path("${meta.id}.seed_permutation_multiqc_jaccard.txt")     , emit: multiqc_jaccard
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    seed_permutation_evaluation.py \\
        --prefix ${meta.id} \\
        --module ${module} \\
        --seeds ${seeds} \\
        --permuted_modules ${permuted_modules} \\
        --permuted_seeds ${permuted_seeds} \\
        --network ${network}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
