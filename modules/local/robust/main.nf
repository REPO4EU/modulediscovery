
process ROBUST {
    label 'process_single'
    container 'djskelton/robust:cc669c6'

    input:
    path seeds
    path network
    // TODO: check whether the values in args should be converted to nextflow values and set defaults via nextflow

    output:
    path "${seeds.baseName}.graphml",        emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''          // Get possible alpha, beta, n, and tau arguments for robust, see TODO above
    """
    python3 /robust/robust.py \\
        $network \\
        $seeds \\
        "${seeds.baseName}.graphml" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
