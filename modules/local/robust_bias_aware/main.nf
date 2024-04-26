process ROBUSTBIASAWARE {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/robust-bias-aware:0.0.1--pyh7cba7a3_1'

    input:
    tuple val(meta), path(seeds)
    val idspace

    output:
    tuple val(meta), path("${seeds.baseName}.graphml")  , emit: module
    path "versions.yml"     , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    """
    robust-bias-aware --seeds $seeds --outfile "${seeds.baseName}.graphml" --namespace $idspace
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        robust-bias-aware: \$(pip show robust-bias-aware | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}
