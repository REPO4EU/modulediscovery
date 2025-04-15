process TOPAS{
    tag "$meta.id"
    label 'process_single'
    container 'community.wave.seqera.io/library/r-dnet_r-base_r-dosnow_r-igraph_pruned:a4e3c1d955d54d45'

    input:
    tuple val(meta), path(seeds), path(network)

    output:
    tuple val(meta), path("${meta.id}.topas.txt")   , emit: module
    path "version.yml"

    when:
    task.ext.when == null || task.ext.when 

    script:
}