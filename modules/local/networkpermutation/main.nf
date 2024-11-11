
process NETWORKPERMUTATION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(network)

    output:
    //tuple val(meta), path("${network.baseName}.*.${network.extension}") , emit: permuted_networks
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    randomize_network.py --network ${network}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
