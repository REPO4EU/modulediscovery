
process PERMUTATIONEVALUATION {
    tag "$meta.id"
    label 'process_single'
    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    tuple val(meta), path(module)
    path(seeds)
    path(permuted_modules)
    path(permuted_seeds)
    path(network)

    output:
    tuple val(meta), path("${meta.id}.permutation_evaluation_summary.tsv")
    tuple val(meta), path("${meta.id}.permutation_evaluation_detailed.tsv")
    tuple val(meta), path("${meta.id}.permutation_multiqc_summary.tsv")     , emit: multiqc_summary
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    permutation_evaluation.py \\
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
