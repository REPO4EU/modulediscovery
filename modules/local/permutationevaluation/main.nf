
process PERMUTATIONEVALUATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(module)
    path(seeds)
    path(permuted_modules)
    path(permuted_seeds)

    output:
    tuple val(meta), path("${meta.id}.permutation_evaluation.tsv")

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    echo "Original module: ${module}
    Original seeds: ${seeds}
    Permuted modules (list): ${permuted_modules}
    Permuted seeds (list): ${permuted_seeds}" > ${meta.id}.permutation_evaluation.tsv
    """
}
