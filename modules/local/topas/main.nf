process TOPAS{
    tag "$meta.id"
    label 'process_single'
    container 'community.wave.seqera.io/library/r-dnet_r-base_r-dosnow_r-igraph_pruned:03e60a99a8a266a6'

    input:
    tuple val(meta), path(seeds), path(network)

    output:
    tuple val(meta), path("${meta.id}.topas.txt")   , emit: module
    path "version.yml"

    when:
    task.ext.when == null || task.ext.when 

    script:
    """
    run_topas.R -n $network -s $seeds -o ${meta.id}.topas.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed s/R //g')
    END_VERSIONS
    """
}