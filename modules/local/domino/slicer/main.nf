
process DOMINO_SLICER {
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/domino:1.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/domino:1.0.0--pyhdfd78af_0' }"

    input:
    path network


    output:
    path "${network.baseName}.slices.txt",      emit: slices
    path "versions.yml",                        emit: versions
    when:
    task.ext.when == null || task.ext.when


    script:
    """
    slicer -n $network -o ${network.baseName}.slices.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
