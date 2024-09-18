process DRUGSTONEEXPORT{
    tag "$meta.id"
    label 'process_single'

    container "docker.io/kerstingjohannes/modulediscovery:1.0.0"

    input:
    tuple val(meta), path(module)
    val(id_space)

    output:
    tuple val(meta), path("${meta.id}.drugstonelink.txt")   , emit: link
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    drugstone.py -m $module -i $id_space -o ${meta.id}.drugstonelink.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
