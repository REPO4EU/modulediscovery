process NETWORKPERMUTATIONEVALUATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(module)
    path(permuted_modules)

    output:
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    echo "
    network_permutation_evaluation.py \\
        --prefix ${meta.id} \\
        --module ${module} \\
        --permuted_modules ${permuted_modules} \\
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
