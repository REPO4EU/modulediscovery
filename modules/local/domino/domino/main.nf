
process DOMINO_DOMINO {     // Process name, should be all upper case. Only the part before "_" will be used to define the output folder
    tag "$meta.id"
    label 'process_low'     // Used to allocate resources, "process_low" uses 2 threads and 12GB memory, for more labels see conf/base.config
    conda "bioconda::domino=1.0.0"  // Define software deployment via conda, "container" is more important right now
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/domino:1.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/domino:1.0.0--pyhdfd78af_0' }"   // The preferred way to two define a container, if a biocontainer is available

    input:
    tuple val(meta), path(seeds), path (network), path(slices)

    output:                                                                       // Define the expected outputs, the "emit:" keyword defines, how the output can be accessed by other processes
    tuple val(meta), path("${seeds.baseName}/modules.out")  , emit: modules       // DOMINO will place a modules.out file in a folder, named as the seed gene file
    path "versions.yml"                                     , emit: versions      // The software versions, in this case it is only the python version


    when:
    task.ext.when == null || task.ext.when

    // Run DOMINO, it supports muliple threads, which are set using the "label" keyword and accessed using $task.cpus
    script:
    def args = task.ext.args ?: ''          // Get possible optional arguments (e.g. turning off visualization, see conf/modules.config)
    """
    domino \\
    $args \\
    -p $task.cpus \\
    -a $seeds \\
    -n $network \\
    -s $slices \\
    -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
