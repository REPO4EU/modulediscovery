process ROBUSTBIASAWARE {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/robust-bias-aware:0.0.1--pyh7cba7a3_1'

    input:
    tuple val(meta), path(seeds), path (network)
    val idspace

    output:
    tuple val(meta), path("${meta.id}.graphml")  , emit: module
    path "versions.yml"     , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''          // Get possible alpha, beta, n, and tau arguments for robust, see TODO above

    """
    robust-bias-aware --seeds $seeds --outfile "${meta.id}.graphml" --namespace $idspace --network $network  $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        robust-bias-aware: \$(pip show robust-bias-aware | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}
