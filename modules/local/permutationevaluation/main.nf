
process PERMUTATIONEVALUATION {
    tag "$meta.id"
    label 'process_single'
    container "community.wave.seqera.io/library/graph-tool_pip_csv.py_l-sys_pruned:e95dbddbf81ad212"

    input:
    tuple val(meta), path(module)
    path(seeds)
    path(permuted_modules)
    path(permuted_seeds)
    path(network)

    output:
    tuple val(meta), path("${meta.id}.permutation_evaluation.tsv")
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    permutation_evaluation.py \\
        --module ${module} \\
        --seeds ${seeds} \\
        --permuted_modules ${permuted_modules} \\
        --permuted_seeds ${permuted_seeds} \\
        --network ${network}

    echo "Original module: ${module}
    Original seeds: ${seeds}
    Permuted modules (list): ${permuted_modules}
    Permuted seeds (list): ${permuted_seeds}" > ${meta.id}.permutation_evaluation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
