process DIAMOND {
    label 'process_single'
    container 'djskelton/diamond:2437974'

    input:
    path seeds
    path network
    val n
    val alpha

    output:
    path "*.txt",        emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python /DIAMOnD/DIAMOnD.py \\
        $network \\
        $seeds \\
        $n \\
        $alpha

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
