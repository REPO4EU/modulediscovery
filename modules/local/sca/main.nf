process SCA{
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(seeds), path (network)

    output:
    tuple val(meta), path("Predicted_genelist.txt")   , emit: module
    path "versions.yml"                               , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    SCA2024.py $network $seeds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
