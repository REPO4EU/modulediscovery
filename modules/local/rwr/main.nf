
// Many additional examples for nextflow modules are available at https://github.com/nf-core/modules/tree/master/modules/nf-core

process RWR {                           // Process name, should be all upper case
    label 'process_single'                  // Used to allocate resources, "process_single" uses one thread and 6GB memory, for labels see conf/base.config
    container 'docker.io/djskelton/diamond:2437974'   // The container on docker hub, other repositories are possible, use conda keyword to set a conda environment

    input:                                  // Define the input channels
    path seeds                              // Path to seeds file
    path network                            // Path to a network file
    val scaling                             // RWR specific parameter "sclaing"
    val symmetrical                         // RWR spefific parameter "symmetrical"
    val r                                   // RWR specific parameter "r"

    output:                                 // Define output files, "emit" is only used to access the corresponding outputs externally
    path "*.txt",        emit: module       // Define a pattern for the output file (can also be the full name, if known), emit -> the active module
    path "versions.yml", emit: versions     // Software versions, this is not essential but nice, the collected versions will be part of the final multiqc report

    when:
    task.ext.when == null || task.ext.when  // Allows to prevent the execution of this process via a workflow logic, just put it in

    // The script for executint RWR, in this case, the .py script is shipped with the container
    // Access inputs, parameters, etc. with the "$" operator
    // The part starting with "cat <<-END_VERSIONS > versions.yml" only collects software versions for the versions.yml file, not essential
    script:
    """
    python ../../../bin/rwr.py \\
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
