
// Many additional examples for nextflow modules are available at https://github.com/nf-core/modules/tree/master/modules/nf-core

process DIAMOND {                           // Process name, should be all upper case
    tag "$meta.id"                          // Used to display the process in the progress overview
    label 'process_single'                  // Used to allocate resources, "process_single" uses one thread and 6GB memory, for labels see conf/base.config
    container 'docker.io/djskelton/diamond:2437974'   // The container on docker hub, other repositories are possible, use conda keyword to set a conda environment

    input:                                  // Define the input channels
    tuple val(meta), path(seeds)            // Path to seeds file
    path network                            // Path to a network file
    val n                                   // DIAMOnD specific parameter "n"
    val alpha                               // DIAMOnD spefific parameter "alpha"

    output:                                 // Define output files, "emit" is only used to access the corresponding outputs externally
    tuple val(meta), path("*.txt"),     emit: module       // Define a pattern for the output file (can also be the full name, if known), emit -> the active module
    path "versions.yml",                emit: versions     // Software versions, this is not essential but nice, the collected versions will be part of the final multiqc report

    when:
    task.ext.when == null || task.ext.when  // Allows to prevent the execution of this process via a workflow logic, just put it in

    // The script for executint DIAMOnD, in this case, the .py script is shipped with the container
    // Access inputs, parameters, etc. with the "$" operator
    // The part starting with "cat <<-END_VERSIONS > versions.yml" only collects software versions for the versions.yml file, not essential
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
