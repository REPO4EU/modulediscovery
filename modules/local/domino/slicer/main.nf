
process DOMINO_SLICER {             // Process name, should be all upper case. Only the part before "_" will be used to define the output folder
    label 'process_single'          // Used to allocate resources, "process_single" uses one thread and 6GB memory, for labels see conf/base.config
    conda "bioconda::domino=1.0.0"  // Define software deployment via conda, "container" is more important right now
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/domino:1.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/domino:1.0.0--pyhdfd78af_0' }" // The preferred way to two define a container, if a biocontainer is available

    input:
    tuple val(meta), path (network)                       // The input network

    output:                                                // Define the expected outputs, the "emit:" keyword defines, how the output can be accessed by other processes
    tuple val(meta), path("${meta.id}.slices.txt"), emit: slices
    path "versions.yml"                           , emit: versions


    when:
    task.ext.when == null || task.ext.when


    // The script for executint slicer
    // Access inputs, parameters, etc. with the "$" operator
    // The part starting with "cat <<-END_VERSIONS > versions.yml" only collects software versions for the versions.yml file, not essential
    script:
    """
    slicer -n $network -o ${meta.id}.slices.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
