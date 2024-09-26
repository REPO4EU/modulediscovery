
process RWR {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    tuple val(meta), path(seeds), path (network)    // Input files
    val scaling                                     // RWR specific parameter "scaling"
    val symmetrical                                 // RWR spefific parameter "symmetrical"
    val r                                           // RWR specific parameter "r"

    output:
    tuple val(meta), path("*.txt") , emit: module
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rwr.py \\
        $network \\
        $seeds \\
        $scaling \\
        $symmetrical \\
        $r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
