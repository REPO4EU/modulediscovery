process SPD {
    tag "$meta.id"  
    label 'process_single'

    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    tuple val(meta), path(subnetwork)
    path network

    output:
    tuple val(meta), path("${meta.id}.spd.gt"), emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    spd_annotation_tool.py -s $subnetwork -n $network -o "${subnetwork.baseName}_spd.gt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
