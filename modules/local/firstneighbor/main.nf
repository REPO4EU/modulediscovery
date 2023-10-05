
// Many additional examples for nextflow modules are available at https://github.com/nf-core/modules/tree/master/modules/nf-core

process FIRSTNEIGHBOR {                       // Process name, should be all upper case
    label 'process_single'                    // Used to allocate resources, "process_single" uses one thread and 6GB memory, for labels see conf/base.config
    conda 'conda-forge::graph-tool=2.58'
    container 'tiagopeixoto/graph-tool:release-2.58'

    input:                                    // Define the input channels
    path seeds                                // Path to seeds file
    path network                              // Path to a network file

    output:                                   // Define output files, "emit" is only used to access the corresponding outputs externally
    path "${seeds.baseName}.gt", emit: module // Define a pattern for the output file (can also be the full name, if known), emit -> the active module
    path "versions.yml", emit: versions       // Software versions, this is not essential but nice, the collected versions will be part of the final multiqc report

    when:
    task.ext.when == null || task.ext.when  // Allows to prevent the execution of this process via a workflow logic, just put it in

    // The script for executing firstneighbor, in this case, the .py script is in the bin folder
    // Access inputs, parameters, etc. with the "$" operator
    // The part starting with "cat <<-END_VERSIONS > versions.yml" only collects software versions for the versions.yml file, not essential
    script:
    """
    firstneighbor_tool.py -n $network -s $seeds -o "${seeds.baseName}.gt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
