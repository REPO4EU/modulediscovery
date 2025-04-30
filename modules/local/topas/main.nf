process TOPAS{
    tag "$meta.id"
    label 'process_single'
    container 'docker.io/motan04/topas:latest'

    input:
    tuple val(meta), path(seeds), path(network)

    output:
    tuple val(meta), path("${meta.id}.topas.txt")   , emit: module
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when 

    script:
    """
    run_topas.R -n $network -s $seeds -o ${meta.id}.topas.txt 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | grep "^R version" | cut -d " " -f1-3)
    END_VERSIONS
    """
}