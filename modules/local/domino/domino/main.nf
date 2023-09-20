
process DOMINO_DOMINO {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/domino:1.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/domino:1.0.0--pyhdfd78af_0' }"

    input:
    path seeds
    path network
    path slices


    output:
    path "*/modules.out",                       emit: modules
    path "versions.yml",                        emit: versions


    when:
    task.ext.when == null || task.ext.when


    script:
    """
    domino -p $task.cpus -a $seeds -n $network -s $slices -v false -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
