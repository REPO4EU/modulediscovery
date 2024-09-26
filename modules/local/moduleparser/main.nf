process MODULEPARSER {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    tuple val(meta), path(module), path(seeds), path(network)
    val tool

    output:
    tuple val(meta), path("${meta.id}.gt")  , emit: module
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    module_parser.py $network -t $tool -l DEBUG -m $module -s $seeds -o ${meta.id}.gt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
